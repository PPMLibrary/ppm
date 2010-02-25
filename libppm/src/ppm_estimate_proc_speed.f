      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_estimate_proc_speed
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine can be used by the user to estimate
      !                 the relative speeds of the processors (will be used in
      !                 subs2proc for load balancing). The estimation is
      !                 done by computing a number of Lennard-Jones PP
      !                 interactions.
      !
      !  Input        :  mintime       (F) Minimum time for which benchmark
      !                                    is required to run on each
      !                                    processor. OPTIONAL. Default is
      !                                    0.1 seconds.
      !                  maxtime       (F) Benchmark stops as soon as
      !                                    slowest processor has run for
      !                                    this time (provided that mintime
      !                                    requirement is met). OPTIONAL.
      !                                    Default is 5 seconds.
      !                  Npart         (I) Initial number of LJ-particles
      !                                    to use. OPTIONAL. Default is
      !                                    1000. Set this to a smaller
      !                                    value when running on a slow
      !                                    processor. Increase is done
      !                                    automatically to meet the time
      !                                    requirements.
      !
      !  Input/output : 
      !
      !  Output       :  proc_speed(:) (F) Relative speeds of all processors
      !                                    from 0 to ppm_nproc-1. The numbers
      !                                    do sum up to 1.
      !                  info          (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_estimate_proc_speed.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2004/08/30 10:53:10  ivos
      !  Added optional arguments for mintime, maxtime and Npart.
      !
      !  Revision 1.6  2004/07/26 07:45:25  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.5  2004/07/20 16:46:01  ivos
      !  Removed unnecessary cpp ifs.
      !
      !  Revision 1.4  2004/07/16 14:46:24  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/05/06 07:43:53  ivos
      !  bugfix: declaration of cutoff2 fixed.
      !
      !  Revision 1.2  2004/04/05 08:52:32  walther
      !  IF (ldone .EQ. .FALSE) replaced by: IF (.NOT.ldone).
      !
      !  Revision 1.1  2004/03/31 10:47:24  ivos
      !  Initial implementation. Tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_estimate_proc_speed_s(proc_speed,info,mintime,maxtime,  &
     &    Npart)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_estimate_proc_speed_d(proc_speed,info,mintime,maxtime,  &
     &    Npart)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
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
      REAL(MK), DIMENSION(:) , POINTER       :: proc_speed
      REAL(MK), OPTIONAL     , INTENT(IN   ) :: mintime,maxtime
      INTEGER , OPTIONAL     , INTENT(IN   ) :: Npart
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0,lmyeps,tim,tim1,cutoff2
      REAL(MK)                         :: xmax,ymax,zmax,r2,r2i,r6i,fs,en
      REAL(MK)                         :: tmin,tmax
      REAL(MK), DIMENSION(3)           :: rx,ri,rj,fi,fj
      REAL(MK), DIMENSION(:), POINTER  :: alltim
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER  :: sendtim
#endif
      REAL(ppm_kind_double)            :: rsum
      INTEGER                          :: i,j,N,iopt
      INTEGER, DIMENSION(1)            :: ldl,ldu
      CHARACTER(LEN=ppm_char)          :: mesg
      LOGICAL                          :: ldone
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_estimate_proc_speed',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_estimate_proc_speed',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (PRESENT(mintime)) THEN
              IF (mintime .LT. 0.0_MK) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_estimate_proc_speed',  &
     &                'mintime must be >= 0.0',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (PRESENT(maxtime)) THEN
              IF (maxtime .LT. 0.0_MK) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_estimate_proc_speed',  &
     &                'maxtime must be >= 0.0',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (PRESENT(Npart)) THEN
              ! 0 itself is NOT permitted since it would bever increase
              ! (multiplicative growth)
              IF (Npart .LE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_estimate_proc_speed',  &
     &                'Npart must be > 0',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Set timing limits. Stop as soon as slowest processor takes more
      !  than tmax seconds, but require every processor to run for at least
      !  tmin seconds for sufficient statistics.
      !-------------------------------------------------------------------------
      tmax = 5.0_MK
      tmin = 0.1_MK
      IF (PRESENT(maxtime)) tmax = maxtime
      IF (PRESENT(mintime)) tmax = mintime

      !-------------------------------------------------------------------------
      !  Set number of particles to start with. The user can override this if
      !  running on a very slow processor.
      !-------------------------------------------------------------------------
      N    = 1000     
      IF (PRESENT(Npart))      N = Npart

      !-------------------------------------------------------------------------
      !  Allocate proc_speed array
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = 0
      ldu(1) = ppm_nproc-1
      CALL ppm_alloc(proc_speed,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_estimate_proc_speed',    &
     &        'relative processor speeds PROC_SPEED',__LINE__,info)
          GOTO 9999
      ENDIF
      proc_speed = 0.0_MK

      !-------------------------------------------------------------------------
      !  If there only is one processor, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .EQ. 1) THEN
          proc_speed(0) = 1.0_MK
          CALL ppm_write(ppm_rank,'ppm_estimate_proc_speed', &
     &           'Only one processor found. Terminating.',info) 
          GOTO 9999
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate memory for all processor times
      !---------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nproc
      CALL ppm_alloc(alltim,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_estimate_proc_speed',    &
     &        'timings from all processors ALLTIM',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Set constants for benchmark calculation
      !-------------------------------------------------------------------------
      xmax = 1.0_MK
      ymax = 1.0_MK
      zmax = 1.0_MK
      cutoff2 = 0.25_MK
      ldone = .FALSE.

      DO WHILE (.NOT. ldone)
          !---------------------------------------------------------------------
          !  Timed LJ computation
          !---------------------------------------------------------------------
          CALL ppm_util_time(tim)
          en = 0.0_MK
          DO i=1,N-1        ! loop over all pairs of particles
             DO j=i+1,N 
                CALL RANDOM_NUMBER(ri)
                CALL RANDOM_NUMBER(rj)
                rx(1) = rj(1) - ri(1)                 ! vector between i and j
                rx(2) = rj(2) - ri(2)
                rx(3) = rj(3) - ri(3)
                rx(1) = rx(1)-xmax*NINT(rx(1)/xmax)     ! nearest image
                rx(2) = rx(2)-ymax*NINT(rx(2)/ymax)  
                rx(3) = rx(3)-zmax*NINT(rx(3)/zmax)  
                r2 = SUM(rx(:)**2)                      ! square of distance
                IF (r2 .LT. cutoff2) THEN 
                   r2i = 1.0_MK/r2
                   r6i = r2i*r2i*r2i
                   ! force factor using Lennard-Jones
                   fs = 48.0_MK*r2i*r6i*(r6i-0.5_MK)
                   fi(1) = fi(1)+fs*rx(1)       ! add force on particle i
                   fi(2) = fi(2)+fs*rx(2)
                   fi(3) = fi(3)+fs*rx(3)
                   fj(1) = fj(1)-fs*rx(1)       ! add opposite force on particle j
                   fj(2) = fj(2)-fs*rx(2)
                   fj(3) = fj(3)-fs*rx(3)
                   ! update potential energy of system
                   en = en+4.0_MK*r6i*(r6i-1.0_MK)
                ENDIF
             ENDDO
          ENDDO
          CALL ppm_util_time(tim1)
          tim = tim1-tim

          alltim = 0.0_MK
#ifdef __MPI
          !---------------------------------------------------------------------
          !  Broadcast timings to all processors
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = ppm_nproc
          CALL ppm_alloc(sendtim,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_estimate_proc_speed',    &
     &            'send buffer for timings SENDTIM',__LINE__,info)
              GOTO 9999
          ENDIF
          sendtim(1:ppm_nproc) = tim
#if   __KIND == __SINGLE_PRECISION
          CALL MPI_Alltoall(sendtim,1,MPI_REAL,alltim,1,MPI_REAL,ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
          CALL MPI_Alltoall(sendtim,1,MPI_DOUBLE_PRECISION,alltim,1,      &
     &                      MPI_DOUBLE_PRECISION,ppm_comm,info)
#endif
          iopt = ppm_param_dealloc
          CALL ppm_alloc(sendtim,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_estimate_proc_speed',     &
     &            'send buffer for timings SENDTIM',__LINE__,info)
          ENDIF
#else
          alltim(1) = tim
#endif
          
          !---------------------------------------------------------------------
          !  Check if test was sufficient
          !---------------------------------------------------------------------
          DO i=1,ppm_nproc
              IF (alltim(i) .GT. tmax) ldone = .TRUE.
          ENDDO
          DO i=1,ppm_nproc
              IF (alltim(i) .LT. tmin) ldone = .FALSE.
          ENDDO

          !---------------------------------------------------------------------
          !  Go again with doubled time if needed
          !---------------------------------------------------------------------
          IF (.NOT.ldone) THEN
              ! CEILING is needed, because it would never grow for small N
              ! if NINT was used !
              N = CEILING(1.414_MK*REAL(N,MK))
          ENDIF

      ENDDO     ! WHILE(.NOT.ldone)

      !-------------------------------------------------------------------------
      !  Convert timings to relative speeds
      !-------------------------------------------------------------------------
      DO i=0,ppm_nproc-1
          proc_speed(i) = 1.0_MK/alltim(i+1)
      ENDDO
      tim = SUM(proc_speed(0:ppm_nproc-1))
      DO i=0,ppm_nproc-1
          proc_speed(i) = proc_speed(i)/tim
      ENDDO

      !-------------------------------------------------------------------------
      !  Check that numbers add up to 1
      !-------------------------------------------------------------------------
      IF (ABS(SUM(proc_speed(0:ppm_nproc-1)) - 1.0_MK) .GT. lmyeps) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_bad_sum,'ppm_estimate_proc_speed', &
     &        'proc_speed must sum up to 1.0',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Diagnostics output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I6,A)') 'Used ',N,' Lennard-Jones particles.'
          CALL ppm_write(ppm_rank,'ppm_estimate_proc_speed',mesg,info)
          CALL ppm_write(ppm_rank,'ppm_estimate_proc_speed', &
                 '---------------------------------------------',info)
          DO i=0,ppm_nproc-1
              WRITE(mesg,'(A,I4,A,F10.7)') 'rank ',i,' has relative speed ', &
     &                    proc_speed(i)
              CALL ppm_write(ppm_rank,'ppm_estimate_proc_speed',mesg,info)
          ENDDO
          CALL ppm_write(ppm_rank,'ppm_estimate_proc_speed', &
                 '---------------------------------------------',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the internal ppm_proc_speed values ????????
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Deallocate local memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(alltim,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_estimate_proc_speed',     &
     &        'timings from all processors ALLTIM',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_estimate_proc_speed',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_estimate_proc_speed_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_estimate_proc_speed_d
#endif
