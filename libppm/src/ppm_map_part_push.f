      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_map_part_push
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine pushes particle data onto the send buffer.
      !
      !  Input        : pdata(:[,:])(O) particle data. Overloaded for
      !                                 single, double, integer, logical,
      !                                 single complex and double complex.
      !                                 Can be either 1d or 2d array.
      !                 lda         (I) the leading dimension of the pdata.
      !                                 In the case of pdata being a 1d
      !                                 array use lda=1.
      !                                 pdata arrays.
      !                 Npart       (I) the number of particles 
      !                 pushpp      (L) OPTIONAL: if true then pdata is assumed
      !                                 to contain the particle positions (xp)
      !
      !  Input/output : info        (I) return status. 0 on success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. The packing could be performed more efficiently.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_push.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.18  2006/10/24 14:20:29  ivos
      !  bugfix: ibuffer was not initialized correctly (Jens fixed it).
      !
      !  Revision 1.17  2006/10/10 20:43:42  walther
      !  Added the optional argument: pushpp.
      !  If pushpp = .TRUE. then the ghost_offset
      !  is added to the psendbuffer.
      !
      !  Revision 1.16  2006/09/04 18:34:51  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.15  2005/02/16 05:06:28  ivos
      !  Bugfix in the unrolled versions.
      !
      !  Revision 1.14  2005/02/16 02:17:09  ivos
      !  Unrolled loops for 1,2,3,4,5 dimensional vectors to allow
      !  vectorization. Only vector mappings with lda.GT.5 do not
      !  vectorize. Scalar mappings do. Tested on NEC SX-5.
      !
      !  Revision 1.13  2004/10/01 16:09:08  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.12  2004/09/28 00:02:29  walther
      !  Bug fix in the (re)allocation of the ppm_sendbuffer(s/d): before
      !  we assumed we could compute the required size, based on Npart, b
      !  but now we compute it exactly based on the loop information.
      !
      !  Revision 1.11  2004/07/26 07:42:46  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.10  2004/07/19 07:41:53  ivos
      !  Overloaded particle push and pop operations for 1d data arrays.
      !  Added new routines to module and needed checks to interface
      !  routine.
      !
      !  Revision 1.9  2004/05/17 15:48:29  oingo
      !  Changed the check for Npart from .LE. 0 to .LT. 0
      !
      !  Revision 1.8  2004/03/05 13:43:24  ivos
      !  bugfix: index errors corrected for COMPLEX version. COMPLEX version
      !  is now tested.
      !
      !  Revision 1.7  2004/02/24 17:12:40  ivos
      !  Added overloaded versions for single complex and double complex 
      !  particle data.
      !
      !  Revision 1.6  2004/02/05 09:01:43  walther
      !  Moved the #include to column 1.
      !
      !  Revision 1.5  2004/01/26 12:36:52  ivos
      !  Changed ppm_buffer2part from a 2D array (npart,topoid) to a 1D array
      !  (Npart), since it does not depend on the topoid (only one mapping at 
      !  the time can be pending since the sendbuffer itself is only 1D as 
      !  well!).
      !  In map_part_push: removed topoid from the argument list in effect as
      !  it is no longer needed (only 1 mapping can be pending).
      !
      !  Revision 1.4  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.3  2004/01/23 11:31:22  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.2  2003/12/10 09:12:02  walther
      !  Updated the header text.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      ! ATTN: DIM is the dimension of the pdata array and not the space
      ! dimension ppm_dim!
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_push_1ds(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_push_1dd(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_1dsc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_1ddc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_push_1di(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_push_1dl(pdata,lda,Npart,info,pushpp)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_push_2ds(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_push_2dd(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_2dsc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_push_2ddc(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_push_2di(pdata,lda,Npart,info,pushpp)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_push_2dl(pdata,lda,Npart,info,pushpp)
#endif 
#endif 
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == 1
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:  ), INTENT(IN   )    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  ), INTENT(IN   )    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: pdata
#else
      REAL(MK), DIMENSION(:  ), INTENT(IN   )    :: pdata
#endif

#elif __DIM == 2
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:), INTENT(IN   )    :: pdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:), INTENT(IN   )    :: pdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | \
      __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
#else
      REAL(MK), DIMENSION(:,:), INTENT(IN   )    :: pdata
#endif
#endif
      INTEGER                 , INTENT(IN   )    :: lda,Npart
      INTEGER                 , INTENT(  OUT)    :: info
      LOGICAL, OPTIONAL       , INTENT(IN   )    :: pushpp
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ipart,ibuffer,icount
      INTEGER               :: iopt,ldb,incr
      REAL(MK)              :: t0
      LOGICAL               :: lpushpp
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_push',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
#if   __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_push',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == 1
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_push',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_push',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  The default is that we do not push the particle positions
      !-------------------------------------------------------------------------
      IF (PRESENT(pushpp)) THEN
         lpushpp = pushpp
      ELSE
         lpushpp = .FALSE.
      ENDIF 

      !-------------------------------------------------------------------------
      !  Increment the buffer set 
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set + 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the buffer dimension and type
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu(1)  = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_push',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_push',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  A complex number is treated as two reals. Cannot change lda
      !  because it is INTENT(IN)
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ldb = 2*lda
#else
      ldb = lda
#endif 

      !-------------------------------------------------------------------------
      !  Store the dimension and type
      !-------------------------------------------------------------------------
      ppm_buffer_dim(ppm_buffer_set)  = ldb
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#elif __KIND == __INTEGER
      ppm_buffer_type(ppm_buffer_set) = ppm_integer
#elif __KIND == __LOGICAL
      ppm_buffer_type(ppm_buffer_set) = ppm_logical
#endif

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist() 
      !-------------------------------------------------------------------------
      ibuffer = ppm_nsendbuffer

      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer 
         !----------------------------------------------------------------------
         incr = 0 
         DO i=1,ppm_nsendlist
            incr = incr + (ppm_psendbuffer(i+1)-ppm_psendbuffer(i))*ldb
         ENDDO
         ldu(1) = ppm_nsendbuffer + incr
         iopt   = ppm_param_alloc_grow_preserve
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_push',     &
     &           'global send buffer PPM_SENDBUFFERD',__LINE__,info)
             GOTO 9999
         ENDIF

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access the particles belonging to the i-th processor in the 
            !  sendlist
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            !  Store the particle data in the buffer
            !-------------------------------------------------------------------
#if   __DIM == 2
            !-------------------------------------------------------------------
            !  Unrolled for lda=1
            !-------------------------------------------------------------------
            IF (lda .EQ. 1) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),    &
     &                ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=2
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 2) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),    &
     &                ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=3
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 3) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),    &
     &                ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=4
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 4) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(4,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(4,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),    &
     &                ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=5
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 5) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                  ppm_sendbufferd(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(4,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = pdata(5,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =         &
     &                       REAL(pdata(5,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(5,ipart)),  &
     &                ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(1,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(2,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(3,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(4,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(4,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) =     &
     &                       REAL(pdata(5,ipart),ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = AIMAG(pdata(5,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbufferd(ibuffer) = REAL(pdata(1,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(2,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(3,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(4,ipart),    &
     &                ppm_kind_double)
                  ibuffer = ibuffer + 1
                  ppm_sendbufferd(ibuffer) = REAL(pdata(5,ipart),    &
     &                ppm_kind_double)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(5,ipart)) THEN
                     ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                  ELSE
                     ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Not unrolled for the rest. Vector length will be lda!!
            !-------------------------------------------------------------------
            ELSE
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  DO k=1,lda
                     ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                     ppm_sendbufferd(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                     ppm_sendbufferd(ibuffer) = pdata(k,ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =         &
     &                          REAL(pdata(k,ipart),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(k,ipart)),  &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =     &
     &                          REAL(pdata(k,ipart),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = AIMAG(pdata(k,ipart))
#elif  __KIND == __INTEGER
                     ppm_sendbufferd(ibuffer) = REAL(pdata(k,ipart),    &
     &                   ppm_kind_double)
#elif  __KIND == __LOGICAL
                     IF (pdata(k,ipart)) THEN
                        ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                     ELSE
                        ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                     ENDIF
#endif
                  ENDDO
               ENDDO
            ENDIF
#elif  __DIM == 1
            !-------------------------------------------------------------------
            !  Scalar version
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the particle id 
               !----------------------------------------------------------------
               ipart = ppm_buffer2part(j)
               ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
               ppm_sendbufferd(ibuffer) = pdata(ipart)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
               ibuffer = ibuffer + 1
               ppm_sendbufferd(ibuffer) = REAL(AIMAG(pdata(ipart)),  &
     &             ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
               ibuffer = ibuffer + 1
               ppm_sendbufferd(ibuffer) = AIMAG(pdata(ipart))
#elif  __KIND == __INTEGER
               ppm_sendbufferd(ibuffer) = REAL(pdata(ipart),ppm_kind_double)
#elif  __KIND == __LOGICAL
               IF (pdata(ipart)) THEN
                  ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
               ELSE
                  ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
               ENDIF
#endif
            ENDDO
#endif
         ENDDO                ! i=1,ppm_nsendlist
      !-------------------------------------------------------------------------
      !  SINGLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      ELSE
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer 
         !----------------------------------------------------------------------
         incr = 0 
         DO i=1,ppm_nsendlist
            incr = incr + (ppm_psendbuffer(i+1)-ppm_psendbuffer(i))*ldb
         ENDDO
         ldu(1) = ppm_nsendbuffer + incr
         iopt   = ppm_param_alloc_grow_preserve
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_push',     &
     &           'global send buffer PPM_SENDBUFFERS',__LINE__,info)
             GOTO 9999
         ENDIF

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access the particles belonging to the i-th processor in the 
            !  sendlist
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            !  Store the particle data in the buffer
            !-------------------------------------------------------------------
#if   __DIM == 2
            !-------------------------------------------------------------------
            !  Unrolled for lda=1
            !-------------------------------------------------------------------
            IF (lda .EQ. 1) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=2
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 2) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=3
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 3) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=4
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 4) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(4,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(4,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Unrolled for lda=5
            !-------------------------------------------------------------------
            ELSEIF (lda .EQ. 5) THEN
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                  ppm_sendbuffers(ibuffer) = pdata(1,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(2,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(3,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(4,ipart)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = pdata(5,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(1,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(2,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(3,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(4,ipart)),  &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(5,ipart)),  &
     &                ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(1,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(2,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(3,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(4,ipart))
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = AIMAG(pdata(5,ipart))
#elif  __KIND == __INTEGER
                  ppm_sendbuffers(ibuffer) = REAL(pdata(1,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(2,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(3,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(4,ipart),   &
     &                ppm_kind_single)
                  ibuffer = ibuffer + 1
                  ppm_sendbuffers(ibuffer) = REAL(pdata(5,ipart),   &
     &                ppm_kind_single)
#elif  __KIND == __LOGICAL
                  IF (pdata(1,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(2,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(3,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(4,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
                  ibuffer = ibuffer + 1
                  IF (pdata(5,ipart)) THEN
                     ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                  ELSE
                     ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                  ENDIF
#endif
               ENDDO
            !-------------------------------------------------------------------
            !  Not unrolled. Vector length will be lda !!
            !-------------------------------------------------------------------
            ELSE
               DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                  !-------------------------------------------------------------
                  !  Get the particle id 
                  !-------------------------------------------------------------
                  ipart = ppm_buffer2part(j)
                  DO k=1,lda
                     ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                     ppm_sendbuffers(ibuffer) = pdata(k,ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(k,ipart)),  &
     &                   ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = AIMAG(pdata(k,ipart))
#elif  __KIND == __INTEGER
                     ppm_sendbuffers(ibuffer) = REAL(pdata(k,ipart),   &
     &                   ppm_kind_single)
#elif  __KIND == __LOGICAL
                     IF (pdata(k,ipart)) THEN
                        ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                     ELSE
                        ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                     ENDIF
#endif
                  ENDDO
               ENDDO
            ENDIF
#elif  __DIM == 1
            !-------------------------------------------------------------------
            !  Scalar version
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the particle id 
               !----------------------------------------------------------------
               ipart = ppm_buffer2part(j)
               ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
               ppm_sendbuffers(ibuffer) = pdata(ipart)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
               ibuffer = ibuffer + 1
               ppm_sendbuffers(ibuffer) = REAL(AIMAG(pdata(ipart)),  &
     &             ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
               ibuffer = ibuffer + 1
               ppm_sendbuffers(ibuffer) = AIMAG(pdata(ipart))
#elif  __KIND == __INTEGER
               ppm_sendbuffers(ibuffer) = REAL(pdata(ipart),ppm_kind_single)
#elif  __KIND == __LOGICAL
               IF (pdata(ipart)) THEN
                  ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
               ELSE
                  ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
               ENDIF
#endif
            ENDDO
#endif
         ENDDO         ! i=1,ppm_nsendlist
      ENDIF         ! ppm_kind

      !-------------------------------------------------------------------------
      !  If we are pushing particle positions (if lpushpp is true) we need to
      !  add the offset to the particles; pushing the particles can only occur
      !  for 2D arrays in single or double precions: xp is USUALLY stored as
      !  xp(1:ppm_dim,1:Npart)
      !-------------------------------------------------------------------------
      IF (lpushpp) THEN
         !----------------------------------------------------------------------
         !  the particle positions are per construction (by a call to ghost_get) 
         !  ALWAYS stored from ibuffer = 1 to nsendbuffer(2) - 1 and the 
         !  ppm_ghost_offset is therefore also stored from 1 - nsendbuffer(2)-1
         !  thus icount = 1 and ibuffer = ppm_nsendlist (which has not yet been
         !  updated (see below)
         !----------------------------------------------------------------------
         icount  = 0 
         ibuffer = ppm_nsendbuffer
         IF (ppm_kind.EQ.ppm_kind_double) THEN
            DO i=1,ppm_nsendlist
               IF (lda.EQ.2) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=2
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = ppm_sendbufferd(ibuffer) &
     &                                        + ppm_ghost_offsetd(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = ppm_sendbufferd(ibuffer) &
     &                                        + ppm_ghost_offsetd(icount) 
                  ENDDO
               ELSEIF (lda.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=3
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id 
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = ppm_sendbufferd(ibuffer) &
     &                                        + ppm_ghost_offsetd(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = ppm_sendbufferd(ibuffer) & 
     &                                        + ppm_ghost_offsetd(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbufferd(ibuffer) = ppm_sendbufferd(ibuffer) & 
     &                                        + ppm_ghost_offsetd(icount) 
                  ENDDO
               ENDIF ! end of lda = 3
            ENDDO ! enddo of nsendlist
         ELSE ! buffers are single precision 
            DO i=1,ppm_nsendlist
               IF (lda.EQ.2) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=2
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = ppm_sendbuffers(ibuffer) &
     &                                        + ppm_ghost_offsets(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = ppm_sendbuffers(ibuffer) &
     &                                        + ppm_ghost_offsets(icount) 
                  ENDDO
               ELSEIF (lda.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for lda=3
                  !-------------------------------------------------------------
                  DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
                     !----------------------------------------------------------
                     !  Get the particle id 
                     !----------------------------------------------------------
                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = ppm_sendbuffers(ibuffer) &
     &                                        + ppm_ghost_offsets(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = ppm_sendbuffers(ibuffer) & 
     &                                        + ppm_ghost_offsets(icount) 

                     ibuffer = ibuffer + 1
                     icount  = icount  + 1
                     ppm_sendbuffers(ibuffer) = ppm_sendbuffers(ibuffer) & 
     &                                        + ppm_ghost_offsets(icount) 
                  ENDDO
               ENDIF ! end of lda = 3
            ENDDO ! enddo of nsendlist
         ENDIF ! end of single/double precision buffers
      ENDIF ! end of lpushpp

      !-------------------------------------------------------------------------
      !  Update the particle buffer count
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_push',t0,info)
      RETURN
#if   __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_push_1ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_push_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_1ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_push_1di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_push_1dl
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_push_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_push_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_2dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_push_2ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_push_2di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_push_2dl
#endif
#endif
