      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_connect_distrib
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine distributes the connections to all
      !                 processors such that every processor has all the
      !                 connections it has at least one particle for
      !
      !  Input        : lda       (I) leading dimension of cd(:,:)
      !                 id(:)     (I) local to global particle number mapping
      !                 Npart     (I) number of particles
      !
      !  Input/output : cd(:,:)   (I) connection data
      !                 Ncon      (I) number of connections
      !
      !  Output       : info      (I) Return status. 0 on success.
      !
      !  Remarks      : This routine directly sends and receives the
      !                 connection data as integers
      !
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_connect_distrib.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/07/21 09:14:13  kotsalie
      !  FRIDAY
      !
      !  Revision 1.4  2004/11/11 15:24:32  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.3  2004/10/01 16:33:35  ivos
      !  cosmetics.
      !
      !  Revision 1.2  2004/10/01 16:09:03  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.1  2004/07/26 12:02:07  ivos
      !  REnamed due to compiler name size limit.
      !
      !  Revision 1.2  2004/07/26 08:54:23  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.1  2004/07/26 08:24:42  ivos
      !  Check-in after module cleanup. All ppm_map_part_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !  Revision 1.9  2004/07/26 07:42:43  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.8  2004/07/23 17:39:03  oingo
      !  Bugfix: The connection data was copied from/to the send buffer with
      !          the wrong length.
      !  Changed id from POINTER to INTENT(IN)
      !
      !  Revision 1.7  2004/07/20 09:00:02  oingo
      !  The id array is now scalar (i.e. dummy leading dimension removed)
      !
      !  Revision 1.6  2004/07/16 14:10:09  ivos
      !  Clarified comment header.
      !
      !  Revision 1.5  2004/07/14 11:51:32  oingo
      !  Rearranged the #ifdef __MPI such that orphan or doubled connections
      !  will be removed even if __MPI is not defined
      !
      !  Revision 1.4  2004/06/10 16:20:00  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.3  2004/06/02 07:55:13  ivos
      !  did not compile on ifort 8.0. inserted missing newline.
      !
      !  Revision 1.2  2004/05/28 11:51:28  walther
      !  Cosmetics.
      !
      !  Revision 1.1  2004/05/11 13:27:27  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_connect_distrib(cd,lda,Ncon,id,Npart,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_invert_list
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
      INTEGER, DIMENSION(:,:), POINTER       :: cd
      INTEGER                , INTENT(IN   ) :: lda
      INTEGER                , INTENT(INOUT) :: Ncon
      INTEGER, DIMENSION(:)  , INTENT(IN   ) :: id
      INTEGER                , INTENT(IN   ) :: Npart
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,k,l,msend,mrecv,iopt,ipos,icon
      REAL(ppm_kind_double) :: t0
      INTEGER               :: tag,cd_size,ncons_local,itarget,isource
      INTEGER               :: hops,ilda
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpstatus
#endif

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_connect_distrib',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Ncon .LT. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,      &
     &            'ppm_map_connect_distrib',   &
     &            'Number of connections must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Get the inverse of the local to global particle id mapping
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(id_temp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'allocating connection array ID_TEMP',__LINE__,info)
          GOTO 9999
      ENDIF

      id_temp(1:Npart) = id(1:Npart)

      CALL ppm_util_invert_list(id_temp,id_inv,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'inverting id list',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt = ppm_param_dealloc
      CALL ppm_alloc(id_temp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'deallocating connection array ID_TEMP',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Prepare the buffer for the local connections
      !-------------------------------------------------------------------------
      cd_size = 10
      ldu(1) = lda
      ldu(2) = cd_size
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(cd_local,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'allocating connection array CD_LOCAL',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the connections that stay on this processor
      !-------------------------------------------------------------------------
      ncons_local = 0
      DO j = 1,Ncon
         DO i = 1,lda
            IF ((cd(i,j) .GE. LBOUND(id_inv,1)) .AND.  &
     &          (cd(i,j) .LE. UBOUND(id_inv,1))) THEN
                !---------------------------------------------------------------
                !  If any local particle is in the current connection, keep it
                !---------------------------------------------------------------
                IF (id_inv(cd(i,j)) .GT. -HUGE(cd(i,j))) THEN
                    ncons_local = ncons_local + 1

                    !-----------------------------------------------------------
                    !  Increase the buffer if necessary
                    !-----------------------------------------------------------
                    IF (ncons_local .GT. cd_size) THEN
                        cd_size = cd_size + 10
                        iopt = ppm_param_alloc_grow_preserve
                        ldu(1) = lda
                        ldu(2) = cd_size
                        CALL ppm_alloc(cd_local,ldu,iopt,info)
                        IF (info .NE. 0) THEN
                            info = ppm_error_fatal
                            CALL ppm_error(ppm_err_alloc,              &
     &                          'ppm_map_connect_distrib',     &
     &                          'allocating connection array CD_LOCAL',&
     &                          __LINE__,info)
                            GOTO 9999
                        ENDIF
                    ENDIF
                   !------------------------------------------------------------
                   !  Store this connection in the array and go to the next one
                   !------------------------------------------------------------
                   DO ilda=1,lda  
                    cd_local(ilda,ncons_local) = cd(ilda,j)
                   ENDDO  
                   EXIT
                ENDIF
            ENDIF
         ENDDO
      ENDDO

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Compute whom to send the connections and who to receive the
      !  connections from
      !-------------------------------------------------------------------------
      itarget = MOD(ppm_nproc + ppm_rank - 1, ppm_nproc)
      isource = MOD(ppm_rank + 1, ppm_nproc)

      !-------------------------------------------------------------------------
      !  Send it to all processors
      !-------------------------------------------------------------------------
      hops = ppm_nproc - 1

      !-------------------------------------------------------------------------
      !  Calculate the size of the data that we are going to send
      !-------------------------------------------------------------------------
      msend = Ncon*lda

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendbuffer
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = msend
      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'allocating send buffer SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill the sendbuffer with all connections we have
      !-------------------------------------------------------------------------
      ipos = 1
      DO j = 1,Ncon
         DO i = 1,lda
            sendbuffer(ipos) = cd(i,j)
            ipos = ipos + 1
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Send the connection data to the next proc and receive them from the
      !  previous
      !-------------------------------------------------------------------------
      DO k = 1,hops
         !----------------------------------------------------------------------
         !  First tell the neighbour how many connection we are going to send.
         !  The other neighbour tells us how many connections we will get
         !----------------------------------------------------------------------
         tag = 100
         CALL MPI_SendRecv(msend,1,MPI_INTEGER,itarget,tag, &
     &                     mrecv,1,MPI_INTEGER,isource,tag, &
     &                     ppm_comm,mpstatus,info)

         !----------------------------------------------------------------------
         !  Allocate memory for the receivebuffer since we now know its size
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu(1) = mrecv
         CALL ppm_alloc(recvbuffer,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &           'allocating receive buffer RECVBUFFER',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Send and receive the actual connection data
         !----------------------------------------------------------------------
         tag = 200
         CALL MPI_SendRecv(sendbuffer,msend,MPI_INTEGER,itarget,tag, &
     &                     recvbuffer,mrecv,MPI_INTEGER,isource,tag, &
     &                     ppm_comm,mpstatus,info)

         !----------------------------------------------------------------------
         !  Pick out the connections that we need
         !----------------------------------------------------------------------
         DO j = 1,mrecv,lda
            DO i = 1,lda
               IF ((recvbuffer(j+(i-1)) .GE. LBOUND(id_inv,1)) .AND.  &
     &             (recvbuffer(j+(i-1)) .LE. UBOUND(id_inv,1))) THEN
                   IF (id_inv(recvbuffer(j+(i-1))) .GT. -HUGE(recvbuffer(j+(i-1)))) THEN
                       ncons_local = ncons_local + 1

                       IF (ncons_local .GT. cd_size) THEN
                           cd_size = cd_size + 10
                           iopt = ppm_param_alloc_grow_preserve
                           ldu(1) = lda
                           ldu(2) = cd_size
                           CALL ppm_alloc(cd_local,ldu,iopt,info)
                           IF (info .NE. 0) THEN
                               info = ppm_error_fatal
                               CALL ppm_error(ppm_err_alloc,             &
     &                            'ppm_map_connect_distrib',     &
     &                            'allocating connection array CD_LOCAL',&
     &                            __LINE__,info)
                               GOTO 9999
                           ENDIF
                       ENDIF
                         cd_local(:,ncons_local) = recvbuffer(j:(j+lda-1))
                       EXIT
                   ENDIF
               ENDIF
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !  The just received connections will be send to the next proc
         !----------------------------------------------------------------------
         msend = mrecv

         !----------------------------------------------------------------------
         !  Swap the sendbuffer and receivebuffer
         !----------------------------------------------------------------------
         tempbuffer => sendbuffer
         sendbuffer => recvbuffer
         recvbuffer => tempbuffer
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'deallocating send buffer SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(recvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'deallocating receive buffer RECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Shrink the array if we allocated more than needed
      !-------------------------------------------------------------------------
      IF ((ncons_local .LT. cd_size) .AND. (ncons_local .NE. 0)) THEN
          iopt = ppm_param_alloc_fit_preserve
          ldu(1) = lda
          ldu(2) = ncons_local
          CALL ppm_alloc(cd_local,ldu,iopt,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &           'allocating connection CD_LOCAL',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(id_inv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'deallocating inverted id list ID_INV',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Mark all connections that we have more than once
      !-------------------------------------------------------------------------
      icon = 0
      DO j = 1,ncons_local-1
         IF (cd_local(1,j) .EQ. -1) CYCLE

         DO k = j+1,ncons_local
            l = 0
            DO i = 1,lda
               IF (cd_local(i,j) .EQ. cd_local(i,k)) THEN
                   l = l + 1
               ENDIF
            ENDDO

            IF (l .EQ. lda) THEN
                cd_local(1,k) = -1
                icon = icon + 1
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Prepare the final connection array
      !-------------------------------------------------------------------------
      ldu(1) = lda
      ldu(2) = ncons_local - icon
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(cd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'allocating connection array CD',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Go through all connections and save all that are not marked
      !-------------------------------------------------------------------------
      Ncon = 0
      DO j = 1,ncons_local
         IF (cd_local(1,j) .NE. -1) THEN
             Ncon = Ncon + 1
             DO ilda=1,lda
              cd(ilda,Ncon) = cd_local(ilda,j)
             ENDDO  
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(cd_local,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_distrib',  &
     &        'deallocating connection array CD_LOCAL',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_connect_distrib',t0,info)
      RETURN

      END SUBROUTINE ppm_map_connect_distrib
