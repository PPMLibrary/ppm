      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_connect_send
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine sends and receives the connection data
      !                 according to the information stored in ppm_buffer2part
      !
      !  Input        : lda       (I) : The leading dimension of cd(:,:)
      !                 id(:)     (I) : local to global particle number mapping
      !                 Npart     (I) : number of particles
      !
      !  Input/output : cd(:,:)   (I) : connection data
      !                 Ncon      (I) : number of connections
      !
      !  Output       : info      (I) : Return status. 0 on success.
      !
      !  Remarks      : This routine directly sends and receives the
      !                 connection data as integers
      !
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_connect_send.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/11/11 15:24:32  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.5  2004/10/01 16:33:35  ivos
      !  cosmetics.
      !
      !  Revision 1.4  2004/10/01 16:09:03  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.3  2004/07/26 11:46:54  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 08:54:24  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.1  2004/07/26 08:24:44  ivos
      !  Check-in after module cleanup. All ppm_map_part_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !  Revision 1.8  2004/07/26 07:42:44  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/07/23 17:39:13  oingo
      !  Changed id from POINTER to INTENT(IN)
      !
      !  Revision 1.6  2004/07/20 08:59:35  oingo
      !  The id array is now scalar (i.e. dummy leading dimension removed)
      !
      !  Revision 1.5  2004/07/19 12:28:05  ivos
      !  bugfix: Added newline before INCLUDE. pgf90 was confused.
      !
      !  Revision 1.4  2004/07/19 11:16:16  oingo
      !  Added check for an empty buffer
      !
      !  Revision 1.3  2004/07/14 11:52:09  oingo
      !  Cosmetics
      !
      !  Revision 1.2  2004/06/10 16:20:01  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/05/11 13:27:46  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_connect_send(cd,lda,Ncon,id,Npart,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_part, ONLY: ppm_target_topoid
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
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,l,msend,mrecv,iopt,ipos,part,qpart,icon
      REAL(ppm_kind_double) :: t0
      INTEGER               :: tag,cd_size,ncons_local
      INTEGER               :: hops
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
#endif

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_connect_send',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Ncon .LT. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,'ppm_map_connect_send', &
     &            'Number of connections must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      IF (ppm_target_topoid .LT. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_nomap,'ppm_map_connect_send',  &
     &        'No mapping defined',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if the send buffer is empty
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set .LT. 1) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_buffer_empt,'ppm_map_connect_send',  &
     &         'No mapping defined!',__LINE__,info)
          GOTO 9999
      ENDIF

#ifndef __MPI
      !-------------------------------------------------------------------------
      !  If we only have one processor there is not much to do
      !-------------------------------------------------------------------------
      GOTO 9999
#else
      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(csend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating array CSEND',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(crecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating array CRECV',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = ppm_nsendlist
      ldu(2) = Ncon
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating array PSEND',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = Npart
      CALL ppm_alloc(id_send,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating list ID_SEND',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the arrays and counter
      !-------------------------------------------------------------------------
      id_send(:) = 0
      csend(:)   = 0
      crecv(:)   = 0
      psend(:,:) = 0

      msend = 0
      mrecv = 0

      !-------------------------------------------------------------------------
      !  Calculate how many particles will stay on this processor
      !-------------------------------------------------------------------------
      ipos = ppm_psendbuffer(2) - ppm_psendbuffer(1)

      !-------------------------------------------------------------------------
      !  Count how many connections we have to send to the other processors
      !-------------------------------------------------------------------------
      DO k = 2,ppm_nsendlist
         qpart = ppm_psendbuffer(k+1) - ppm_psendbuffer(k)
         !----------------------------------------------------------------------
         !  Loop over all particles that go the k-th processor
         !----------------------------------------------------------------------
         DO i = ipos+1,(ipos+qpart)
            part = id(ppm_buffer2part(i))
            id_send(ppm_buffer2part(i)) = part
            !-------------------------------------------------------------------
            !  Loop over all connections to find the particle part
            !-------------------------------------------------------------------
            DO j = 1,Ncon
               DO l = 1,lda
                  IF (part .EQ. cd(l,j)) THEN
                      !---------------------------------------------------------
                      !  Count and mark this particle if it is not already
                      !  counted and marked (this particles goes to the k-th
                      !  processor)
                      !---------------------------------------------------------
                      IF (psend(k,j) .EQ. 0) THEN
                          csend(k) = csend(k) + 1
                          psend(k,j) = 1
                          !----------------------------------------------------
                          !  Mark all particles that will be send away. This
                          !  is needed to do a faster check for orphan
                          !  connections
                          !----------------------------------------------------
                          psend(1,j) = 1
                      ENDIF
                      EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !  Communicate the number of connections we have to send and to
         !  receive
         !----------------------------------------------------------------------
         tag = 100
         CALL MPI_SendRecv(csend(k),1,MPI_INTEGER,ppm_isendlist(k),tag, &
     &                     crecv(k),1,MPI_INTEGER,ppm_irecvlist(k),tag, &
     &                     ppm_comm,status,info)

         msend = msend + csend(k)
         mrecv = mrecv + crecv(k)

         ipos = ipos + qpart
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the connection data array in order to hold the incoming
      !  connections
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = lda
      ldu(2) = Ncon + mrecv
      CALL ppm_alloc(cd_local,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating connectin array CD_LOCAL',__LINE__,info)
          GOTO 9999
      ENDIF

      cd_local(:,:) = 0
      cd_local(:,1:Ncon) = cd(:,:)

      !-------------------------------------------------------------------------
      !  Allocate memory for the send buffer and receive buffer
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = MAXVAL(csend,1)*lda
      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating send buffer SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = MAXVAL(crecv,1)*lda
      CALL ppm_alloc(recvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'allocating receive buffer RECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Send and receive the connections
      !-------------------------------------------------------------------------
      icon = 1
      DO k = 2,ppm_nsendlist
         !----------------------------------------------------------------------
         !  Fill the sendbuffer for the k-th processor
         !----------------------------------------------------------------------
         ipos = 1
         DO j = 1,Ncon
            IF (psend(k,j) .EQ. 1) THEN
                DO i = 1,lda
                   sendbuffer(ipos) = cd_local(i,j)
                   ipos = ipos + 1
                ENDDO
            ENDIF
         ENDDO

         tag = 200
         CALL MPI_SendRecv(sendbuffer,csend(k)*lda,MPI_INTEGER,      &
     &                     ppm_isendlist(k),tag,                     &
     &                     recvbuffer,crecv(k)*lda,MPI_INTEGER,      &
     &                     ppm_irecvlist(k),tag,ppm_comm,status,info)

         !----------------------------------------------------------------------
         !  Add the received connections to the connection data array
         !----------------------------------------------------------------------
         IF (crecv(k) .NE. 0) THEN
             ipos = 1
             DO j = 1,crecv(k)
                DO i = 1,lda
                   cd_local(i,Ncon+icon) = recvbuffer(ipos)
                   ipos = ipos + 1
                ENDDO
                icon = icon + 1
             ENDDO
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  We now have a new number of connections
      !-------------------------------------------------------------------------
      ncons_local = Ncon + mrecv

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(csend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating array CSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(crecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating array CRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating send buffer SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(recvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating receive buffer RECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Invert list of particles that has been sent away
      !-------------------------------------------------------------------------
      CALL ppm_util_invert_list(id_send,id_inv,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'inverting id list',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt = ppm_param_dealloc
      CALL ppm_alloc(id_send,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating temporary id list ID_SEND',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find and mark all orphan connections
      !-------------------------------------------------------------------------
      icon = 0
      DO j = 1,Ncon
         IF (psend(1,j) .EQ. 1) THEN
             l = 0
             DO i = 1,lda
                IF ((cd_local(i,j) .GE. LBOUND(id_inv,1)) .AND.  &
     &              (cd_local(i,j) .LE. UBOUND(id_inv,1))) THEN
                    IF (id_inv(cd_local(i,j)) .GT. -HUGE(cd_local(i,j))) THEN
                        l = l + 1
                    ENDIF
                ELSE
                    l = l + 1
                ENDIF
             ENDDO
             IF (l .EQ. lda) THEN
                 cd_local(1,j) = -1
                 icon = icon + 1
             ENDIF
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Mark all connections that we have more than once
      !-------------------------------------------------------------------------
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
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
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
             cd(:,Ncon) = cd_local(:,j)
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating array PSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(cd_local,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_connect_send',  &
     &        'deallocating connection array CD_LOCAL',__LINE__,info)
          GOTO 9999
      ENDIF

#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_connect_send',t0,info)
      RETURN

      END SUBROUTINE ppm_map_connect_send
