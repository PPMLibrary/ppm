      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_io_open
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine opens an IO channel for parallel
      !                 reading and/or writing and associates the given
      !                 I/O unit with the file indicated.
      !
      !  Input        : filename       (C) name of the file to which IO
      !                                    refers
      !                 actn           (I) IO action. One of:
      !                                       ppm_param_io_write
      !                                       ppm_param_io_read
      !                                       ppm_param_io_read_write
      !                 posn           (I) Position for IO in file. This is
      !                                    only relavent for writing
      !                                    action. One of:
      !                                       ppm_param_io_replace
      !                                       ppm_param_io_append
      !                 pfmt           (I) Format of IO. One of:
      !                                       ppm_param_io_ascii
      !                                       ppm_param_io_binary
      !                 pmode          (I) parallel mode. Specifies whether
      !                                    every processor writes its data
      !                                    to a local file or if all data
      !                                    is routed to rank 0, which is
      !                                    writing everything into one
      !                                    single file. One of:
      !                                       ppm_param_io_distributed
      !                                       ppm_param_io_centralized
      !
      !  Input/output : iUnit          (I) Fortran IO unit to be used. If
      !                                    -1 is given on input, the
      !                                    routine will find an unused unit
      !                                    and return it here.
      !
      !  Output       : info           (I) return code. 0 on success.
      !
      !  Remarks      : There are lots of MPI_Allreduce and MPI_Bcast going
      !                 on for the processors to negotiate a unit number
      !                 and catch all errors. There must be a cheaper way
      !                 of doing this...
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_open.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2004/10/01 16:33:34  ivos
      !  cosmetics.
      !
      !  Revision 1.8  2004/10/01 16:09:01  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.7  2004/07/26 15:38:47  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.6  2004/07/26 07:45:26  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.5  2004/07/16 14:46:26  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.4  2004/06/10 13:45:23  ivos
      !  bugfix: inserted missing initialization of lopen.
      !
      !  Revision 1.3  2004/05/14 14:45:26  ivos
      !  Tested. Added argument check for iUnit .LE. 0 and Allreduce of info
      !  after the actual OPEN (since this can fail, too).
      !
      !  Revision 1.2  2004/05/13 11:41:19  ivos
      !  bugfix: all processors check units now in the distributed case and
      !  only if the unit is available on all of them, it is taken.
      !
      !  Revision 1.1  2004/05/06 10:44:58  ivos
      !  Preliminary check-in. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_open(iUnit,filename,actn,posn,pfmt,pmode,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_alloc
      USE ppm_module_io_unused_unit
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
       INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(IN   ) :: actn,posn,pfmt,pmode
      CHARACTER(LEN=*)       , INTENT(IN   ) :: filename
      INTEGER                , INTENT(INOUT) :: iUnit
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lopen,loall
      CHARACTER(LEN=ppm_char)          :: mesg
      CHARACTER(LEN=12)                :: cactn,cfmt,cposn,cstat
      INTEGER                          :: maxused,iopt,info2
      INTEGER, DIMENSION(1)            :: ldc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_open',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_open',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Decode and check arguments
      !-------------------------------------------------------------------------
      IF (iUnit .LE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_open',    &
     &        'Unit number needs to be > 0',__LINE__,info)
          GOTO 9999
      ENDIF

      IF     (actn .EQ. ppm_param_io_read) THEN
          cactn = 'READ'
      ELSEIF (actn .EQ. ppm_param_io_write) THEN
          cactn = 'WRITE'
      ELSEIF (actn .EQ. ppm_param_io_read_write) THEN
          cactn = 'READWRITE'
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_io_open',    &
     &       'Invalid action type specified.',__LINE__,info)
         GOTO 9999
      ENDIF

      IF     (pfmt .EQ. ppm_param_io_ascii) THEN
          cfmt = 'FORMATTED'
      ELSEIF (pfmt .EQ. ppm_param_io_binary) THEN
          cfmt = 'UNFORMATTED'
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_io_open',    &
     &       'Invalid format type specified.',__LINE__,info)
         GOTO 9999
      ENDIF

      IF     ((pmode .NE. ppm_param_io_distributed) .AND.    &
     &        (pmode .NE. ppm_param_io_centralized)) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_io_open',    &
     &       'Invalid I/O mode specified.',__LINE__,info)
         GOTO 9999
      ENDIF

      IF (actn .EQ. ppm_param_io_read) THEN
          cstat = 'OLD'
          cposn = 'REWIND'
      ELSE
          IF     (posn .EQ. ppm_param_io_append) THEN
              cstat = 'UNKNOWN'
              cposn = 'APPEND'
          ELSEIF (posn .EQ. ppm_param_io_replace) THEN
              cstat = 'REPLACE'
              cposn = 'REWIND'
          ELSE
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_io_open',    &
     &           'Invalid position type specified.',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that file is not already connected to another unit
      !-------------------------------------------------------------------------
      lopen = .FALSE.
      IF (pmode .EQ. ppm_param_io_centralized) THEN
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(FILE=filename,OPENED=lopen)
          ENDIF
#ifdef __MPI
          !---------------------------------------------------------------------
          !  Allreduce of lopen
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(lopen,loall,1,MPI_LOGICAL,MPI_LOR,   &
     &            ppm_comm,info)
          lopen = loall
#endif
          IF (lopen) THEN
             info = ppm_error_error
             WRITE(mesg,'(3A)') 'File ',TRIM(filename),' is already opened!'
             CALL ppm_error(ppm_err_open,'ppm_io_open',mesg,  &
     &           __LINE__,info)
             iUnit = -1
             GOTO 9999
          ENDIF
      ELSE
          INQUIRE(FILE=filename,OPENED=lopen)
          IF (lopen) THEN
             info = ppm_error_error
             WRITE(mesg,'(3A)') 'File ',TRIM(filename),' is already opened!'
             CALL ppm_error(ppm_err_open,'ppm_io_open',mesg,  &
     &           __LINE__,info)
             iUnit = -1
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  CENTRALIZED I/O
      !-------------------------------------------------------------------------
      IF (pmode .EQ. ppm_param_io_centralized) THEN
          !---------------------------------------------------------------------
          !  Root finds the next free unit
          !---------------------------------------------------------------------
          IF (iUnit .LE. 0) THEN
              iUnit = 9
              lopen = .TRUE.
              DO WHILE (lopen)
                  IF (ppm_rank .EQ. 0) THEN
                      CALL ppm_io_unused_unit(iUnit,info)
                  ENDIF
#ifdef __MPI
                  !-------------------------------------------------------------
                  !  Bcast to all
                  !-------------------------------------------------------------
                  CALL MPI_Bcast(iUnit,1,MPI_INTEGER,0,ppm_comm,info)
#endif
                  !-------------------------------------------------------------
                  !  If no unused unit was found, all processors fail
                  !-------------------------------------------------------------
                  IF (iUnit .LT. 0) THEN
                     info = ppm_error_error
                     CALL ppm_error(ppm_err_open,'ppm_io_open',    &
     &                   'No I/O unit available.',__LINE__,info)
                     GOTO 9999
                  ENDIF

                  !-------------------------------------------------------------
                  !  All check that no such ppm unit is open locally
                  !-------------------------------------------------------------
                  lopen = .FALSE.
                  IF (ASSOCIATED(ppm_io_unit)) THEN
                      IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
                          IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
                      ENDIF
                  ENDIF
#ifdef __MPI
                  !-------------------------------------------------------------
                  !  Allreduce of lopen
                  !-------------------------------------------------------------
                  CALL MPI_Allreduce(lopen,loall,1,MPI_LOGICAL,MPI_LOR,   &
     &                ppm_comm,info)
                  lopen = loall
#endif
              ENDDO  

          !---------------------------------------------------------------------
          !  Everyone checks the ppm unit, Root the F90 unit in addition
          !---------------------------------------------------------------------
          ELSE
              lopen = .FALSE.
              IF (ppm_rank .EQ. 0) THEN
                  INQUIRE(iUnit,OPENED=lopen)
                  IF (iUnit .EQ. ppm_stdout) lopen = .TRUE.
                  IF (iUnit .EQ. ppm_stderr) lopen = .TRUE.
                  IF (iUnit .EQ. ppm_logfile) lopen = .TRUE.
              ENDIF
              IF (ASSOCIATED(ppm_io_unit)) THEN
                  IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
                      IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
                  ENDIF
              ENDIF
#ifdef __MPI
              !-----------------------------------------------------------------
              !  Allreduce of lopen
              !-----------------------------------------------------------------
              CALL MPI_Allreduce(lopen,loall,1,MPI_LOGICAL,MPI_LOR,   &
     &                ppm_comm,info)
              lopen = loall
#endif
              IF (lopen) THEN
                 info = ppm_error_error
                 WRITE(mesg,'(A,I4,A)') 'Requested I/O unit ',iUnit,  &
     &               ' is already in use. Close other instance first!'
                 CALL ppm_error(ppm_err_unit_open,'ppm_io_open',mesg,  &
     &               __LINE__,info)
                 iUnit = -1
                 GOTO 9999
              ENDIF
          ENDIF
      !-------------------------------------------------------------------------
      !  DISTRIBUTED I/O
      !-------------------------------------------------------------------------
      ELSE
          IF (iUnit .LE. 0) THEN
              iUnit = 9
              lopen = .TRUE.
              DO WHILE (lopen)
                  CALL ppm_io_unused_unit(iUnit,info)
                  !-------------------------------------------------------------
                  !  Fail if no unused unit was found
                  !-------------------------------------------------------------
                  IF (iUnit .LT. 0) THEN
                     info = ppm_error_error
                     CALL ppm_error(ppm_err_open,'ppm_io_open',    &
     &                   'No I/O unit available.',__LINE__,info)
                     GOTO 9999
                  ENDIF
                  !-------------------------------------------------------------
                  !  All check that no such ppm unit is open locally
                  !-------------------------------------------------------------
                  lopen = .FALSE.
                  IF (ASSOCIATED(ppm_io_unit)) THEN
                      IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
                          IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
                      ENDIF
                  ENDIF
              ENDDO  
          ELSE
              lopen = .FALSE.
              INQUIRE(iUnit,OPENED=lopen)
              IF (iUnit .EQ. ppm_stdout) lopen = .TRUE.
              IF (iUnit .EQ. ppm_stderr) lopen = .TRUE.
              IF (iUnit .EQ. ppm_logfile) lopen = .TRUE.
              IF (ASSOCIATED(ppm_io_unit)) THEN
                  IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
                      IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
                  ENDIF
              ENDIF
              IF (lopen) THEN
                 info = ppm_error_error
                 WRITE(mesg,'(A,I4,A)') 'Requested I/O unit ',iUnit,  &
     &               ' is already in use. Close other instance first!'
                 CALL ppm_error(ppm_err_unit_open,'ppm_io_open',mesg,  &
     &               __LINE__,info)
                 iUnit = -1
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Open the file
      !-------------------------------------------------------------------------
      IF (pmode .EQ. ppm_param_io_centralized) THEN
          !---------------------------------------------------------------------
          !  Only root opens the file if centralized IO is requested
          !---------------------------------------------------------------------
          IF (ppm_rank .EQ. 0) THEN
              OPEN(iUnit,FILE=filename,FORM=cfmt,POSITION=cposn,STATUS=cstat, &
     &            ACTION=cactn,IOSTAT=info)
          ENDIF
#ifdef __MPI
          !-----------------------------------------------------------------
          !  Allreduce of info
          !-----------------------------------------------------------------
          CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,   &
     &        ppm_comm,info)
          info = info2
#endif
      ELSE
          !---------------------------------------------------------------------
          !  Otherwise every processor opens a file (of different name
          !  maybe)
          !---------------------------------------------------------------------
          OPEN(iUnit,FILE=filename,FORM=cfmt,POSITION=cposn,STATUS=cstat, &
     &           ACTION=cactn,IOSTAT=info)
      ENDIF
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_open',    &
     &       'File does not exist or no write permissions.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Unit number book-keeping
      !-------------------------------------------------------------------------
      IF (.NOT. ASSOCIATED(ppm_io_unit)) THEN
          iopt = ppm_param_alloc_fit
          ldc(1) = iUnit
          CALL ppm_alloc(ppm_io_unit,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_open',    &
     &            'I/O unit list PPM_IO_UNIT',__LINE__,info)
              GOTO 9999
          ENDIF
          ! all units are unused
          ppm_io_unit(1:iUnit) = 0
      ELSE
          maxused = SIZE(ppm_io_unit,1)
          iopt = ppm_param_alloc_grow_preserve
          ldc(1) = iUnit
          CALL ppm_alloc(ppm_io_unit,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_open',    &
     &            'I/O unit list PPM_IO_UNIT',__LINE__,info)
              GOTO 9999
          ENDIF
          iopt = SIZE(ppm_io_unit,1)
          ! set new units to ununsed
          IF (iopt .GT. maxused) THEN
              ppm_io_unit(maxused+1:iopt) = 0
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Grow lists
      !-------------------------------------------------------------------------
      maxused = MAXVAL(ppm_io_unit)
      maxused = maxused + 1
      iopt = ppm_param_alloc_grow_preserve
      ldc(1) = maxused
      CALL ppm_alloc(ppm_io_mode,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_io_open',    &
     &        'I/O mode list PPM_IO_MODE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_io_format,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_io_open',    &
     &        'I/O format list PPM_IO_FORMAT',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store internal unit number, mode and format
      !-------------------------------------------------------------------------
      ppm_io_unit(iUnit)     = maxused
      ppm_io_mode(maxused)   = pmode
      ppm_io_format(maxused) = pfmt

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I4,A)') 'I/O unit ',iUnit,' opened.'
          CALL ppm_write(ppm_rank,'ppm_io_open',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_open',t0,info)
      RETURN
      END SUBROUTINE ppm_io_open
