      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_doio
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs the actual I/O. It is never
      !                 called by the user program directly, but via
      !                 ppm_io. It only handles 1d arrays. Everything else
      !                 is converted to a 1d array by the calling routine
      !                 ppm_io.
      !
      !  Input        : iUnit          (I) I/O unit (as returned by
      !                                    ppm_io_open).
      !                 actn           (I) I/O action. One of:
      !                                       ppm_param_io_read
      !                                       ppm_param_io_write
      !                 dist           (I) I/O distribution. 
      !                                    For read one of:
      !                                       ppm_param_io_same
      !                                       ppm_param_io_split
      !                                    For write one of:
      !                                       ppm_param_io_root
      !                                       ppm_param_io_concat
      !                                       ppm_param_io_sum
      !                 iofmt          (C) Format string (for ASCII I/O).
      !                                    List-directed I/O is used if
      !                                    string of length 0 is given.
      !                 ioprec         (I) Precision for binary I/O.
      !                                    One of:
      !                                       ppm_param_io_single
      !                                       ppm_param_io_double
      !                                       ppm_param_undefined
      !                                    If undefined is given, no
      !                                    conversion is done.
      !                 lchar          (L) .TRUE. if the integer vector
      !                                    actually corresponds to a
      !                                    character string. In this case
      !                                    the string will be converted
      !                                    back to characters before it is
      !                                    written out. .FALSE. if it
      !                                    is a true INTEGER vector.
      !                 mptag          (I) MPI tag to be used for this I/O.
      !
      !  Input/output : adata(:)       (O) data to be written/read.
      !                                    Overloaded types are: single,
      !                                    double, integer, single complex,
      !                                    double complex, logical,
      !                                    character. Must be a 1d array.
      !
      !  Output       : info           (I) return code. 0 on success.
      !
      !  Remarks      : Uses non-blocking MPI to overlap MPI-IO and
      !                 file-IO. For this reason it needs two buffers. If
      !                 memory is a problem, using only one and blocking
      !                 MPI might be better.
      !
      !                 This routine uses MPI_Irecv and thus, all buffer
      !                 arrays need to be ALLOCATABLEs and not POINTERs!
      !                 Using Fortran90 POINTER structures in conjunction
      !                 with non-blocking MPI_Irecv does NOT WORK! (see
      !                 http://parallel.ru/docs/Parallel/mpi2/node236.html)
      !
      !                 This routine makes use of MPI_Cancel to cancel
      !                 pending send and receive actions when an I/O error
      !                 occured on Root. Canceling a send action is
      !                 expensive and only supported in versions newer than
      !                 mpich 1.2.0.
      !
      !  References   : http://parallel.ru/docs/Parallel/mpi2/node236.html
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_doio.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.15  2006/09/28 21:53:07  ivos
      !  Replaced all the array notation with DO loops. Reason: array copy
      !  is done on the stack, thus crashing the program when large data
      !  are in/output. DO loops have no limits and hence we can now also
      !  read and write BIG restart files...
      !
      !  Revision 1.14  2005/06/12 01:01:51  ivos
      !  bugfix: mdata was not initialized in the distributed IO case.
      !
      !  Revision 1.13  2004/11/11 15:24:30  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.12  2004/10/01 16:08:57  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.11  2004/09/01 10:48:39  ivos
      !  Change: The warning that is issued if not all the data in a record is
      !  actually read does no longer lead to abort on all processors.
      !
      !  Revision 1.10  2004/07/26 15:38:46  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.9  2004/07/26 07:45:24  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.8  2004/06/11 10:21:58  ivos
      !  Changed reclen from I8B=SELECTED_INT_KIND(18) to regular INTEGER
      !  since I8B caused problems on the NEC.
      !
      !  Revision 1.7  2004/06/03 16:06:47  ivos
      !  bugfix: data was not passed out in READ SPLIT case on a single
      !  processor. fixed.
      !
      !  Revision 1.6  2004/06/03 07:52:24  ivos
      !  bugfix: READ SPLIT did not work on only one processor.
      !
      !  Revision 1.5  2004/06/02 13:40:51  ivos
      !  bugfix: label 500 removed from inside IF branch to its beginning.
      !  NEC SX5 compiler complained.
      !
      !  Revision 1.4  2004/06/02 07:13:50  ivos
      !  bugfix: serial version did not compile due to misplaced ifdef in
      !  variable declarations.
      !
      !  Revision 1.3  2004/05/27 15:26:41  ivos
      !  Bugfix: Record completion for READ SAME and READ SPLIT was checked on
      !  all processors. Should only be done by root to avoid unnecessary
      !  warnings.
      !
      !  Revision 1.2  2004/05/27 14:47:40  ivos
      !  READ SAME was extended to support cross-record reading. The routine
      !  was tested.
      !
      !  Revision 1.1  2004/05/27 10:42:50  ivos
      !  Initial implementation. This contains the actual I/O logic and was
      !  extracted from the old version of ppm_io. The atomic read and write
      !  actions have been moved into subroutines.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_doio_s(iUnit,adata,actn,dist,iofmt,ioprec,lchar,  &
     &    mptag,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_doio_d(iUnit,adata,actn,dist,iofmt,ioprec,lchar,  &
     &    mptag,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_doio_sc(iUnit,adata,actn,dist,iofmt,ioprec,lchar, & 
     &    mptag,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_doio_dc(iUnit,adata,actn,dist,iofmt,ioprec,lchar, &
     &    mptag,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_doio_i(iUnit,adata,actn,dist,iofmt,ioprec,lchar,  &
     &    mptag,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_doio_l(iUnit,adata,actn,dist,iofmt,ioprec,lchar,  &
     &    mptag,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_io_read_ascii
      USE ppm_module_io_read_binary
      USE ppm_module_io_write_ascii
      USE ppm_module_io_write_binary
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
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
      INTEGER                          , INTENT(IN   ) :: iUnit
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(INOUT) :: adata
#endif
      INTEGER                          , INTENT(  OUT) :: info
      INTEGER                          , INTENT(IN   ) :: actn,dist,ioprec,mptag
      CHARACTER(LEN=*)                 , INTENT(IN   ) :: iofmt
      LOGICAL                          , INTENT(IN   ) :: lchar
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: reclen
      REAL(MK)                         :: t0
      LOGICAL                          :: lopen
      INTEGER                          :: ifmt,imode
      INTEGER                          :: info2,info3
      CHARACTER(LEN=ppm_char)          :: mesg
      INTEGER, DIMENSION(5)            :: ldc,lda
      INTEGER, DIMENSION(2)            :: ldl
      INTEGER                          :: i,iopt,mdata,ndata,nrec
      INTEGER                          :: ibuffer,jbuffer,ii
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)       , POINTER :: rbuf
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)       , POINTER :: rbuf
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)       , POINTER :: rbuf
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)       , POINTER :: rbuf
#endif
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpstat
      INTEGER                          :: MPTYPE,iproc,rem
      INTEGER                          :: mpreq,jproc
      INTEGER, DIMENSION(:  ), POINTER :: asize
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)       , POINTER :: abuffer,bbuffer
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)       , POINTER :: abuffer,bbuffer
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)       , POINTER :: abuffer,bbuffer
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)       , POINTER :: abuffer,bbuffer
#endif
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_doio',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_doio',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the size of the data
      !-------------------------------------------------------------------------
      mdata = SIZE(adata,1)

      !-------------------------------------------------------------------------
      !  Check if the unit number is open
      !-------------------------------------------------------------------------
      lopen = .FALSE.
      IF (ASSOCIATED(ppm_io_unit)) THEN
          IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
              IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
          ENDIF
      ENDIF
      IF (.NOT. lopen) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_no_unit,'ppm_doio',    &
     &          'Open unit first using ppm_io_open.',__LINE__,info)
      ENDIF
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  If one fails, all have to fail
      !-------------------------------------------------------------------------
      CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm,info3)
      IF (info3 .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &        'MPI_ALLREDUCE of info.',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &          'Unit not open on all processors.',__LINE__,info)
          GOTO 9999
      ENDIF
#else
      IF (info .NE. 0) GOTO 9999
#endif

      !-------------------------------------------------------------------------
      !  Get format and mode for the unit
      !-------------------------------------------------------------------------
      ifmt = ppm_io_format(ppm_io_unit(iUnit))
      imode = ppm_io_mode(ppm_io_unit(iUnit))

#ifndef __MPI
      !-------------------------------------------------------------------------
      !  For the serial version, all IO is distributed.
      !-------------------------------------------------------------------------
      imode = ppm_param_io_distributed
#endif

      !-------------------------------------------------------------------------
      !  Set pointer to read buffer
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION 
      rbuf => rbuf_s
#elif __KIND == __DOUBLE_PRECISION
      rbuf => rbuf_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX 
      rbuf => rbuf_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      rbuf => rbuf_dc
#elif __KIND == __INTEGER
      rbuf => rbuf_i
#elif __KIND == __LOGICAL
      rbuf => rbuf_l
#endif

      !-------------------------------------------------------------------------
      !  DISTRIBUTED OR ONLY ROOT
      !-------------------------------------------------------------------------
      IF ((imode .EQ. ppm_param_io_distributed) .OR.              &
     &    ((dist .EQ. ppm_param_io_root) .AND. (ppm_rank .EQ. 0))) THEN
          INQUIRE(iUnit,OPENED=lopen)
          IF (.NOT. lopen) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_unit,'ppm_doio',    &
     &            'No open file associated with unit.',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (actn .EQ. ppm_param_io_read) THEN
              !-----------------------------------------------------------------
              !  Read from file
              !-----------------------------------------------------------------
              IF (ifmt .EQ. ppm_param_io_binary) THEN
                  ibuffer = 0 
                  DO WHILE (ibuffer .LT. mdata)
                      !-----------------------------------------------------
                      !  Read record length from file
                      !-----------------------------------------------------
                      READ(iUnit,ERR=500,IOSTAT=info) reclen
                      IF (info .NE. 0) GOTO 500
                      nrec = reclen
                      !-----------------------------------------------------
                      !  Allocate record buffer
                      !-----------------------------------------------------
                      iopt = ppm_param_alloc_fit
                      ldl(1) = nrec
                      CALL ppm_alloc(rbuf,ldl,iopt,info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_alloc,'ppm_doio',  &
     &                        'record buffer RBUF',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      jbuffer = 0
                      !-----------------------------------------------------
                      !  Read next record into rbuf
                      !-----------------------------------------------------
                      CALL ppm_io_read_binary(iUnit,rbuf,nrec,ioprec,   &
     &                    lchar,info)
                      IF (info .NE. 0) GOTO 500
                      !-----------------------------------------------------
                      !  Fill bbuffer and update counters
                      !-----------------------------------------------------
                      IF (ibuffer+nrec .LE. mdata) THEN
                          ! adata can hold whole record
                          DO ii=1,nrec
                              adata(ibuffer+ii) = rbuf(ii)
                          ENDDO
                          ibuffer = ibuffer + nrec
                          jbuffer = nrec
                      ELSE
                          ! adata too small. keep remainder
                          jbuffer = ndata - ibuffer
                          DO ii=1,jbuffer
                              adata(ibuffer+ii) = rbuf(ii)
                          ENDDO
                          ibuffer = ndata
                      ENDIF
                  ENDDO      ! WHILE ibuffer.LT.ndata
                  !---------------------------------------------------------
                  !  Check if all data was sold
                  !---------------------------------------------------------
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      IF (jbuffer .NE. nrec) THEN
                          info3 = info
                          info = ppm_error_warning
                          CALL ppm_error(ppm_err_data_miss,'ppm_doio',  &
     &                       'Last record was not completely read',     &
     &                       __LINE__,info)
                          ! restore old info as warnings need not be
                          ! broadcasted
                          info = info3
                      ENDIF
                  ENDIF
              ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                  CALL ppm_io_read_ascii(iUnit,adata,iofmt,lchar,info)
              ENDIF
          ELSEIF (actn .EQ. ppm_param_io_write) THEN
              !-----------------------------------------------------------------
              !  Write data to file
              !-----------------------------------------------------------------
              IF (ifmt .EQ. ppm_param_io_binary) THEN
                  CALL ppm_io_write_binary(iUnit,adata,ioprec,lchar,info)
              ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                  CALL ppm_io_write_ascii(iUnit,adata,iofmt,lchar,info)
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Check if action was successful
          !---------------------------------------------------------------------
 500      IF (info .NE. 0) THEN
              info = ppm_error_error
              GOTO 9999
          ENDIF
          
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  CENTRALIZED
      !-------------------------------------------------------------------------
      ELSEIF ((imode .EQ. ppm_param_io_centralized) .AND.   &
     &        (dist .NE. ppm_param_io_root)) THEN

          !---------------------------------------------------------------------
          !  Set buffer pointer
          !---------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION 
          abuffer => abuffer_s
          bbuffer => bbuffer_s
#elif __KIND == __DOUBLE_PRECISION
          abuffer => abuffer_d
          bbuffer => bbuffer_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX 
          abuffer => abuffer_sc
          bbuffer => bbuffer_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
          abuffer => abuffer_dc
          bbuffer => bbuffer_dc
#elif __KIND == __INTEGER
          abuffer => abuffer_i
          bbuffer => bbuffer_i
#elif __KIND == __LOGICAL
          abuffer => abuffer_l
          bbuffer => bbuffer_l
#endif

          !---------------------------------------------------------------------
          !  Barrier to ensure that output from different time steps is not
          !  mixed
          !---------------------------------------------------------------------
          CALL MPI_Barrier(ppm_comm,info)

          !---------------------------------------------------------------------
          !  Root checks that the file is open
          !---------------------------------------------------------------------
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(iUnit,OPENED=lopen)
              IF (.NOT. lopen) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_no_unit,'ppm_doio',    &
     &                'No open file associated with unit.',__LINE__,info)
              ENDIF
          ENDIF
          !---------------------------------------------------------------------
          !  If one fails, all have to fail
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm,info3)
          IF (ppm_debug .GT. 0) THEN
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'Rank 0 has no access to physical file.',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          !---------------------------------------------------------------------
          !  Determine MPI data type
          !---------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
          MPTYPE = MPI_REAL
#elif __KIND == __SINGLE_PRECISION_COMPLEX
          MPTYPE = MPI_COMPLEX
#elif __KIND == __DOUBLE_PRECISION 
          MPTYPE = MPI_DOUBLE_PRECISION
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
          MPTYPE = MPI_DOUBLE_COMPLEX
#elif __KIND == __INTEGER
          MPTYPE = MPI_INTEGER
#elif __KIND == __LOGICAL
          MPTYPE = MPI_LOGICAL
#endif
          !---------------------------------------------------------------------
          !  Get the data sizes from all other processors
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldl(1) = 0
          ldc(1) = ppm_nproc-1
          CALL ppm_alloc(asize,ldl,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_doio',    &
     &            'processor data size list ASIZE',__LINE__,info)
              GOTO 9999
          ENDIF

          CALL MPI_Gather(mdata,1,MPI_INTEGER,asize,1,MPI_INTEGER, &
     &        0,ppm_comm,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &            'MPI_GATHER of data sizes',__LINE__,info)
              GOTO 9999
          ENDIF
          !---------------------------------------------------------------------
          !  Debug output
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 1 .AND. ppm_rank .EQ. 0) THEN
              WRITE(mesg,'(A)') '  Rank  |             Data size   '
              CALL ppm_write(ppm_rank,'ppm_doio',mesg,info)
              WRITE(mesg,'(A)') '--------------------------------------------'
              CALL ppm_write(ppm_rank,'ppm_doio',mesg,info)
              DO i=0,ppm_nproc-1
                  WRITE(mesg,'(A,I4,A,I18)') '  ',i,'  |   ',asize(i)
                  CALL ppm_write(ppm_rank,'ppm_doio',mesg,info)
              ENDDO
              WRITE(mesg,'(A)') '--------------------------------------------'
              CALL ppm_write(ppm_rank,'ppm_doio',mesg,info)
          ENDIF

          !---------------------------------------------------------------------
          !  WRITE SUM
          !---------------------------------------------------------------------
          IF (dist .EQ. ppm_param_io_sum) THEN
              !-----------------------------------------------------------------
              !  All data needs to be of the same size if this is a write
              !  action with distribution ppm_param_io_sum
              !-----------------------------------------------------------------
              IF (ppm_rank .EQ. 0) THEN
                  DO iproc=1,ppm_nproc-1
                      IF (asize(iproc) .NE. mdata) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_io_data,'ppm_doio',    &
     &               'All data must be of same size for ppm_param_io_sum', &
     &                       __LINE__,info)
                      ENDIF
                  ENDDO
              ENDIF      
              !-----------------------------------------------------------------
              !  If one fails, all have to fail
              !-----------------------------------------------------------------
              CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm,  &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'Non-matching data found.',__LINE__,info)
                  GOTO 9999
              ENDIF
                  
              !-----------------------------------------------------------------
              !  Allocate receive buffer for the sum
              !-----------------------------------------------------------------
              ALLOCATE(abuffer(mdata),STAT=info2)
              IF (info2 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_doio',    &
     &                'receive buffer ABUFFER',__LINE__,info)
                  GOTO 9999
              ENDIF

              !-----------------------------------------------------------------
              !  Sum data from all processors into receive buffer
              !-----------------------------------------------------------------
              CALL MPI_Reduce(adata,abuffer,mdata,MPTYPE,MPI_SUM,0,ppm_comm, &
     &            info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_REDUCE of data',__LINE__,info)
                  GOTO 9999
              ENDIF

              !-----------------------------------------------------------------
              !  Root writes the summed data
              !-----------------------------------------------------------------
              IF (ppm_rank .EQ. 0) THEN
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      CALL ppm_io_write_binary(iUnit,abuffer,ioprec,lchar,info)
                  ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                      CALL ppm_io_write_ascii(iUnit,abuffer,iofmt,lchar,info)
                  ENDIF
              ENDIF    

              !-----------------------------------------------------------------
              !  If one failed, all have to fail
              !-----------------------------------------------------------------
              CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm, &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'I/O failed on Rank 0.',__LINE__,info)
                  GOTO 9999
              ENDIF

          !---------------------------------------------------------------------
          !  WRITE CONCAT
          !---------------------------------------------------------------------
          ELSEIF (dist .EQ. ppm_param_io_concat) THEN
              IF (ppm_rank .EQ. 0) THEN
                  !-------------------------------------------------------------
                  !  Allocate buffer for local data
                  !-------------------------------------------------------------
                  ALLOCATE(abuffer(mdata),STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_doio',     &
     &                    'buffer ABUFFER',__LINE__,info)
                      GOTO 9999
                  ENDIF
                  abuffer = adata
                  !-------------------------------------------------------------
                  !  Receive data from other processors
                  !-------------------------------------------------------------
                  DO iproc=1,ppm_nproc-1
                      ndata = asize(iproc)
                      !---------------------------------------------------------
                      !  Allocate receive buffer
                      !---------------------------------------------------------
                      ALLOCATE(bbuffer(ndata),STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_alloc,'ppm_doio',     &
     &                        'buffer BBUFFER',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      !---------------------------------------------------------
                      !  Start receiving data from iproc
                      !---------------------------------------------------------
                      CALL MPI_Irecv(bbuffer,ndata,MPTYPE,iproc,mptag,  &
     &                    ppm_comm,mpreq,info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                        'MPI_IRECV of data',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      !---------------------------------------------------------
                      !  Write the data of iproc-1
                      !---------------------------------------------------------
                      IF (ifmt .EQ. ppm_param_io_binary) THEN
                          CALL ppm_io_write_binary(iUnit,abuffer,ioprec,lchar,&
     &                        info)
                      ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                          CALL ppm_io_write_ascii(iUnit,abuffer,iofmt,lchar,  &
     &                        info)
                      ENDIF
                      IF (info .NE. 0) GOTO 700
                      !---------------------------------------------------------
                      !  Wait for MPI receive to complete
                      !---------------------------------------------------------
                      CALL MPI_Wait(mpreq,mpstat,info)
                      !---------------------------------------------------------
                      !  Reallocate and copy buffers
                      !---------------------------------------------------------
                      DEALLOCATE(abuffer,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &                        'buffer ABUFFER',__LINE__,info)
                      ENDIF
                      ALLOCATE(abuffer(ndata),STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_alloc,'ppm_doio',     &
     &                        'buffer ABUFFER',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      abuffer = bbuffer
                      DEALLOCATE(bbuffer,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &                        'buffer BBUFFER',__LINE__,info)
                      ENDIF
                  ENDDO        ! iproc

                  !-------------------------------------------------------------
                  !  Write data of rank nproc-1
                  !-------------------------------------------------------------
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      CALL ppm_io_write_binary(iUnit,abuffer,ioprec,lchar,info)
                  ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                      CALL ppm_io_write_ascii(iUnit,abuffer,iofmt,lchar,info)
                  ENDIF
                  IF (info .NE. 0) GOTO 700
              ELSE        ! rank.NE.root
                  !-------------------------------------------------------------
                  !  Send data to root using non-blocking MPI
                  !-------------------------------------------------------------
                  CALL MPI_Isend(adata,mdata,MPTYPE,0,mptag,ppm_comm,mpreq,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                    'MPI_ISEND of data',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDIF     ! rank
              !-----------------------------------------------------------------
              !  If one failed, all have to fail
              !-----------------------------------------------------------------
 700          CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm, &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'I/O failed on Rank 0.',__LINE__,info)
                  ! Cancel pending send since Root failed and will never
                  ! receive it. This operation is expensive as several MPI
                  ! messages are generated. MPICH supports this only in
                  ! version 1.2.0 and newer.
                  IF (mpreq .NE. MPI_REQUEST_NULL) CALL MPI_Cancel(mpreq,info2)
                  GOTO 9999
              ELSE
                  IF (ppm_rank .NE. 0) CALL MPI_Wait(mpreq,mpstat,info)
              ENDIF

          !---------------------------------------------------------------------
          !  READ SAME
          !---------------------------------------------------------------------
          ELSEIF (dist .EQ. ppm_param_io_same) THEN
              !-----------------------------------------------------------------
              !  All data needs to be at least of the size of the data on
              !  root.
              !-----------------------------------------------------------------
              IF (ppm_rank .EQ. 0) THEN
                  DO iproc=1,ppm_nproc-1
                      IF (asize(iproc) .LT. mdata) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_io_data,'ppm_doio',    &
     &                       'some data buffers too small ppm_param_io_same', &
     &                       __LINE__,info)
                      ENDIF
                  ENDDO
              ENDIF    ! rank .EQ. 0
              !-----------------------------------------------------------------
              !  If one fails, all have to fail
              !-----------------------------------------------------------------
              CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm,  &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'Non-matching data found.',__LINE__,info)
                  GOTO 9999
              ENDIF

              !-----------------------------------------------------------------
              !  Root reads the data from the file
              !-----------------------------------------------------------------
              IF (ppm_rank .EQ. 0) THEN
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      ibuffer = 0 
                      DO WHILE (ibuffer .LT. mdata)
                          !-----------------------------------------------------
                          !  Read record length from file
                          !-----------------------------------------------------
                          READ(iUnit,ERR=600,IOSTAT=info) reclen
                          IF (info .NE. 0) GOTO 600
                          nrec = reclen
                          !-----------------------------------------------------
                          !  Allocate record buffer
                          !-----------------------------------------------------
                          iopt = ppm_param_alloc_fit
                          ldl(1) = nrec
                          CALL ppm_alloc(rbuf,ldl,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_doio',  &
     &                            'record buffer RBUF',__LINE__,info)
                              GOTO 9999
                          ENDIF
                          jbuffer = 0
                          !-----------------------------------------------------
                          !  Read next record into rbuf
                          !-----------------------------------------------------
                          CALL ppm_io_read_binary(iUnit,rbuf,nrec,ioprec,   &
     &                        lchar,info)
                          IF (info .NE. 0) GOTO 600
                          !-----------------------------------------------------
                          !  Fill bbuffer and update counters
                          !-----------------------------------------------------
                          IF (ibuffer+nrec .LE. mdata) THEN
                              ! adata can hold whole record
                              DO ii=1,nrec
                                  adata(ibuffer+ii) = rbuf(ii)
                              ENDDO
                              ibuffer = ibuffer + nrec
                              jbuffer = nrec
                          ELSE
                              ! adata too small. keep remainder
                              jbuffer = ndata - ibuffer
                              DO ii=1,jbuffer
                                  adata(ibuffer+ii) = rbuf(ii)
                              ENDDO
                              ibuffer = ndata
                          ENDIF
                      ENDDO      ! WHILE ibuffer.LT.ndata
                      !---------------------------------------------------------
                      !  Check if all data was sold
                      !---------------------------------------------------------
                      IF (ifmt .EQ. ppm_param_io_binary) THEN
                          IF (jbuffer .NE. nrec) THEN
                              info3 = info
                              info = ppm_error_warning
                              CALL ppm_error(ppm_err_data_miss,'ppm_doio',  &
     &                           'Last record was not completely read',     &
     &                           __LINE__,info)
                              ! restore old info as warnings need not be
                              ! broadcasted
                              info = info3
                          ENDIF
                      ENDIF
                  ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                      CALL ppm_io_read_ascii(iUnit,adata,iofmt,lchar,info)
                  ENDIF
              ENDIF    ! rank.EQ.0
              !-----------------------------------------------------------------
              !  If one failed, all have to fail
              !-----------------------------------------------------------------
 600          CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm, &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'I/O failed on Rank 0.',__LINE__,info)
                  GOTO 9999
              ENDIF
              !-----------------------------------------------------------------
              !  Data is broadcast to all others
              !-----------------------------------------------------------------
              CALL MPI_Bcast(adata,mdata,MPTYPE,0,ppm_comm,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_BCAST of data',__LINE__,info)
                  GOTO 9999
              ENDIF

          !---------------------------------------------------------------------
          !  READ SPLIT
          !---------------------------------------------------------------------
          ELSEIF (dist .EQ. ppm_param_io_split) THEN
              IF (ppm_rank .EQ. 0) THEN
                  jproc = MIN(ppm_nproc-1,1)
                  ndata = asize(jproc)
                  !-------------------------------------------------------------
                  !  Allocate read buffer
                  !-------------------------------------------------------------
                  ALLOCATE(bbuffer(ndata),STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_doio',     &
     &                    'buffer BBUFFER',__LINE__,info)
                      GOTO 9999
                  ENDIF
                  !-------------------------------------------------------------
                  !  Read first data chunk 
                  !-------------------------------------------------------------
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      ibuffer = 0 
                      DO WHILE (ibuffer .LT. ndata)
                          !-----------------------------------------------------
                          !  Read record length from file
                          !-----------------------------------------------------
                          READ(iUnit,ERR=900,IOSTAT=info) reclen
                          IF (info .NE. 0) GOTO 900
                          nrec = reclen
                          !-----------------------------------------------------
                          !  Allocate record buffer
                          !-----------------------------------------------------
                          iopt = ppm_param_alloc_fit
                          ldl(1) = nrec
                          CALL ppm_alloc(rbuf,ldl,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_doio',  &
     &                            'record buffer RBUF',__LINE__,info)
                              GOTO 9999
                          ENDIF
                          jbuffer = 0
                          !-----------------------------------------------------
                          !  Read next record into rbuf
                          !-----------------------------------------------------
                          CALL ppm_io_read_binary(iUnit,rbuf,nrec,ioprec,   &
     &                        lchar,info)
                          IF (info .NE. 0) GOTO 900
                          !-----------------------------------------------------
                          !  Fill bbuffer and update counters
                          !-----------------------------------------------------
                          IF (ibuffer+nrec .LE. ndata) THEN
                              ! abuffer can hold whole record
                              DO ii=1,nrec
                                  bbuffer(ibuffer+ii) = rbuf(ii)
                              ENDDO
                              ibuffer = ibuffer + nrec
                              jbuffer = nrec
                          ELSE
                              ! abuffer too small. keep remainder
                              jbuffer = ndata - ibuffer
                              DO ii=1,jbuffer
                                  bbuffer(ibuffer+ii) = rbuf(ii)
                              ENDDO
                              ibuffer = ndata
                          ENDIF
                      ENDDO      ! WHILE ibuffer.LT.ndata
                  ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                      CALL ppm_io_read_ascii(iUnit,bbuffer,iofmt,lchar,info)
                      IF (info .NE. 0) GOTO 900
                  ENDIF
                  !-------------------------------------------------------------
                  !  If we only have one processor, this is it.
                  !-------------------------------------------------------------
                  IF (ppm_nproc .EQ. 1) adata = bbuffer
                  !-------------------------------------------------------------
                  !  Read data for the other processors
                  !-------------------------------------------------------------
                  DO iproc=1,ppm_nproc-1
                      !---------------------------------------------------------
                      !  Processor for which we do file I/O
                      !---------------------------------------------------------
                      jproc = iproc + 1
                      IF (jproc .EQ. ppm_nproc) jproc = 0
                      !---------------------------------------------------------
                      !  Start sending data to iproc
                      !---------------------------------------------------------
                      CALL MPI_Isend(bbuffer,ndata,MPTYPE,iproc,mptag,  &
     &                    ppm_comm,mpreq,info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                        'MPI_ISEND of data',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      !---------------------------------------------------------
                      !  (Re-)allocate read buffer for jproc
                      !---------------------------------------------------------
                      ndata = asize(jproc)
                      ALLOCATE(abuffer(ndata),STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_fatal
                          CALL ppm_error(ppm_err_alloc,'ppm_doio',    &
     &                        'receive buffer ABUFFER',__LINE__,info)
                          GOTO 9999
                      ENDIF
                      !---------------------------------------------------------
                      !  Read data for jproc from file
                      !---------------------------------------------------------
                      IF (ifmt .EQ. ppm_param_io_binary) THEN
                          !-----------------------------------------------------
                          !  Fill in what remains from the last record
                          !-----------------------------------------------------
                          ibuffer = 0
                          IF (jbuffer .LT. nrec) THEN
                              rem = nrec - jbuffer
                              IF (rem .LE. ndata) THEN
                                  ! abuffer large enough to hold all of it
                                  DO ii=1,rem
                                      abuffer(ii) = rbuf(jbuffer+ii)
                                  ENDDO
                                  ibuffer = rem
                                  jbuffer = nrec
                              ELSE
                                  ! abuffer too small. keep remainder
                                  DO ii=1,ndata
                                      abuffer(ii) = rbuf(jbuffer+ii)
                                  ENDDO
                                  jbuffer = jbuffer + ndata
                                  ibuffer = ndata
                              ENDIF
                          ENDIF

                          !-----------------------------------------------------
                          !  Read record until enough data to fill abuffer
                          !  has been read
                          !-----------------------------------------------------
                          DO WHILE (ibuffer .LT. ndata)
                              !-------------------------------------------------
                              !  Read record length from file
                              !-------------------------------------------------
                              READ(iUnit,ERR=900,IOSTAT=info) reclen
                              IF (info .NE. 0) GOTO 900
                              nrec = reclen
                              !-------------------------------------------------
                              !  Allocate record buffer
                              !-------------------------------------------------
                              iopt = ppm_param_alloc_fit
                              ldl(1) = nrec
                              CALL ppm_alloc(rbuf,ldl,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,'ppm_doio',  &
     &                                'record buffer RBUF',__LINE__,info)
                                  GOTO 9999
                              ENDIF
                              jbuffer = 0
                              !-------------------------------------------------
                              !  Read next record into rbuf
                              !-------------------------------------------------
                              CALL ppm_io_read_binary(iUnit,rbuf,nrec,ioprec, &
     &                            lchar,info)
                              IF (info .NE. 0) GOTO 900
                              !-------------------------------------------------
                              !  Fill abuffer and update counters
                              !-------------------------------------------------
                              IF (ibuffer+nrec .LE. ndata) THEN
                                  ! abuffer can hold whole record
                                  DO ii=1,nrec
                                      abuffer(ibuffer+ii) = rbuf(ii)
                                  ENDDO
                                  ibuffer = ibuffer + nrec
                                  jbuffer = nrec
                              ELSE
                                  ! abuffer too small. keep remainder
                                  jbuffer = ndata - ibuffer
                                  DO ii=1,jbuffer
                                      abuffer(ibuffer+ii) = rbuf(ii)
                                  ENDDO
                                  ibuffer = ndata
                              ENDIF
                          ENDDO      ! WHILE ibuffer.LT.ndata
                      ELSEIF (ifmt .EQ. ppm_param_io_ascii) THEN
                          CALL ppm_io_read_ascii(iUnit,abuffer,iofmt,lchar,info)
                          IF (info .NE. 0) GOTO 900
                      ENDIF
                      !---------------------------------------------------------
                      !  Wait for MPI send to complete
                      !---------------------------------------------------------
                      CALL MPI_Wait(mpreq,mpstat,info)
                      !---------------------------------------------------------
                      !  Reallocate and copy buffers
                      !---------------------------------------------------------
                      DEALLOCATE(bbuffer,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &                        'buffer BBUFFER',__LINE__,info)
                      ENDIF
                      IF (jproc .NE. 0) THEN
                          ALLOCATE(bbuffer(ndata),STAT=info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_doio',     &
     &                            'buffer BBUFFER',__LINE__,info)
                              GOTO 9999
                          ENDIF
                      ENDIF
                      IF (jproc .NE. 0) THEN
                          bbuffer = abuffer
                      ELSE
                          adata = abuffer
                      ENDIF
                      DEALLOCATE(abuffer,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &                        'buffer ABUFFER',__LINE__,info)
                      ENDIF
                  ENDDO        ! iproc

                  !-------------------------------------------------------------
                  !  Check if all data was sold
                  !-------------------------------------------------------------
                  IF (ifmt .EQ. ppm_param_io_binary) THEN
                      IF (jbuffer .NE. nrec) THEN
                          info3 = info
                          info = ppm_error_warning
                          CALL ppm_error(ppm_err_data_miss,'ppm_doio',  &
     &                      'Last record was not completely read',__LINE__,info)
                          ! restore old info as warnings need not be
                          ! broadcasted
                          info = info3
                      ENDIF
                  ENDIF
              ELSE        ! rank.NE.root
                  !-------------------------------------------------------------
                  !  Receive data from root using non-blocking MPI
                  !-------------------------------------------------------------
                  CALL MPI_Irecv(adata,mdata,MPTYPE,0,mptag,ppm_comm,mpreq,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                    'MPI_IRECV of data',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDIF     ! rank

              !-----------------------------------------------------------------
              !  If one failed, all have to fail
              !-----------------------------------------------------------------
 900          CALL MPI_Allreduce(info,info2,1,MPI_INTEGER,MPI_SUM,ppm_comm, &
     &            info3)
              IF (info3 .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',    &
     &                'MPI_ALLREDUCE of info.',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (info2 .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mpi_term,'ppm_doio',    &
     &                'I/O failed on Rank 0.',__LINE__,info)
                  ! Cancel pending recv since Root failed and will never
                  ! send anything. This cancel is unproblematic
                  IF (mpreq .NE. MPI_REQUEST_NULL) CALL MPI_Cancel(mpreq,info2)
                  GOTO 9999
              ELSE
                  IF (ppm_rank .NE. 0) CALL MPI_Wait(mpreq,mpstat,info)
              ENDIF

          ENDIF    ! dist
   
          !---------------------------------------------------------------------
          !  Deallocate data size list
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(asize,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &            'processor data size list ASIZE',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Barrier to ensure that output from different time steps is not
          !  mixed
          !---------------------------------------------------------------------
          CALL MPI_Barrier(ppm_comm,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_mpi_fail,'ppm_doio',  &
     &            'MPI barrier failed',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF              ! imode

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      !-------------------------------------------------------------------------
      !  Deallocate buffers
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(sbuffer,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_doio',    &
     &         'single buffer SBUFFER',__LINE__,info)
      ENDIF
      CALL ppm_alloc(dbuffer,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_doio',    &
     &         'double buffer DBUFFER',__LINE__,info)
      ENDIF
#ifdef __MPI
      IF (ppm_rank .EQ. 0) THEN
          iopt = ppm_param_dealloc
          CALL ppm_alloc(rbuf,ldc,iopt,info2)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_doio',     &
     &            'record buffer RBUF',__LINE__,info)
          ENDIF
      ENDIF
      IF (ASSOCIATED(abuffer)) DEALLOCATE(abuffer,STAT=info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_doio',    &
     &         'buffer ABUFFER',__LINE__,info)
      ENDIF
      IF (ASSOCIATED(bbuffer)) DEALLOCATE(bbuffer,STAT=info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_doio',    &
     &         'buffer BBUFFER',__LINE__,info)
      ENDIF
#endif
      CALL substop('ppm_doio',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_doio_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_doio_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_doio_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_doio_dc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_doio_i
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_doio_l
#endif
