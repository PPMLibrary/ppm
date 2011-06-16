      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_copy_image_subs
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine copies the sub domains to ghost sub domains
      !                 in the case of periodic boundary condition and stores 
      !                 the original ID of the sub domain - this to allow the
      !                 find neighbours to search (without wrapping space) the
      !                 list of sub domains and still retrieve the real ID of 
      !                 sub domain. 
      !
      !  Input        : min_phys(:)  (F) : the min. extent of the physical domain
      !                 max_phys(:)  (F) : the amx. extent of the physical domain
      !                 bcdef(:)     (I) : boundary condition definition
      !                 nsubs        (I) : the total number of sub domains
      !
      !  Input/output : min_sub(:,:) (F) : the min. extent of the sub domain 
      !                                    including after the call to this 
      !                                    routine the ghost sub domains
      !                 max_sub(:,:) (F) : the max. extent of the sub domain
      !                                    including after the call to this 
      !                                    routine the ghost sub domains
      !                 subid(:)     (I) : the real id of the sub domain    
      !
      !  Output       : nsubsplus    (I) : the total number of sub domains plus
      !                                    the ghost sub domains
      !                 info         (I) : return status
      !
      !  Remarks      : This is a happy going routine: we are comparing the 
      !                 equality of floating point numbers BUT it is WORKS 
      !                 since the value of subs are ASSIGNED from the phys 
      !                 min/max and NOT computed (byte copy).
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_copy_image_subs.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2006/04/20 15:30:32  ivos
      !  Buffer size is now increased multiplicatively. Better scaling.
      !
      !  Revision 1.2  2006/03/28 19:37:35  ivos
      !  Removed unnecessary alloc.
      !
      !  Revision 1.1  2004/07/26 08:55:31  ivos
      !  Renamed.
      !
      !  Revision 1.10  2004/07/26 07:42:37  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.9  2004/07/06 13:31:21  ivos
      !  Indentation cosmetics.
      !
      !  Revision 1.8  2004/04/02 10:14:06  walther
      !  Bugfix: replaced wrong bcdef(1:3) with bcdef(1,3,5) and added check for
      !  consistency of periodic boundary conditions.
      !
      !  Revision 1.7  2004/01/27 12:15:58  ivos
      !  Bugfix: Added IF(ppm_dim.GT.2) around argument check for third dimension.
      !
      !  Revision 1.6  2004/01/26 14:08:47  walther
      !  Update the header again (min/max_sub is modified by the routine).
      !
      !  Revision 1.5  2004/01/26 13:57:46  walther
      !  Updated the header; check of input arguments; error handling after each
      !  ppm_alloc() call.
      !
      !  Revision 1.4  2004/01/23 17:24:14  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.3  2003/12/12 15:53:15  ivos
      !  Updated comment header.
      !
      !  Revision 1.2  2003/12/10 17:05:11  walther
      !  Added comments in header.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_copy_image_subs_s(min_phys,max_phys,bcdef, &
     &   min_sub,max_sub,nsubs,subid,nsubsplus,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_copy_image_subs_d(min_phys,max_phys,bcdef, &
     &   min_sub,max_sub,nsubs,subid,nsubsplus,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: bcdef
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER , DIMENSION(:)  , POINTER       :: subid
      INTEGER                 , INTENT(IN   ) :: nsubs
      INTEGER                 , INTENT(  OUT) :: nsubsplus,info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      REAL(MK):: t0
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER :: i,j,k
      INTEGER :: istat,iopt,isize
      CHARACTER(ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_copy_image_subs',t0,info)

      !-------------------------------------------------------------------------
      !  Check the input arguments (if in debugging mode)
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         !----------------------------------------------------------------------
         !  The physical boundary is sane ?
         !----------------------------------------------------------------------
          IF (max_phys(1) .LE. min_phys(1)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &            'max_phys(1) must be > min_phys(1)',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (max_phys(2) .LE. min_phys(2)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &            'max_phys(2) must be > min_phys(2)',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_dim .GT. 2) THEN
              IF (max_phys(3) .LE. min_phys(3)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',&
     &                'max_phys(3) must be > min_phys(3)',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF

         !----------------------------------------------------------------------
         !  The periodic BC are sane ?
         !----------------------------------------------------------------------
         IF (bcdef(1).EQ.ppm_param_bcdef_periodic.AND. &
     &       bcdef(2).NE.ppm_param_bcdef_periodic) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'periodic on only one x face !',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (bcdef(3).EQ.ppm_param_bcdef_periodic.AND. &
     &       bcdef(4).NE.ppm_param_bcdef_periodic) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'periodic on only one y face !',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (ppm_dim.GT.2) THEN 
            IF (bcdef(5).EQ.ppm_param_bcdef_periodic.AND. &
     &          bcdef(6).NE.ppm_param_bcdef_periodic) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                        'periodic on only one z face !',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the extend of the physical system 
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         len_phys(k) = max_phys(k) - min_phys(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the current size of subid
      !-------------------------------------------------------------------------
      isize = SIZE(subid)

      !-------------------------------------------------------------------------
      !  Set constants
      !-------------------------------------------------------------------------
      nsubsplus = nsubs
      iopt      = ppm_param_alloc_fit_preserve

      !-------------------------------------------------------------------------
      !  Add ghost domains for the periodic system
      !-------------------------------------------------------------------------
      IF (bcdef(1).EQ.ppm_param_bcdef_periodic) THEN
         !----------------------------------------------------------------------
         !  Copy in the first direction
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            !-------------------------------------------------------------------
            !  Two dimensions
            !-------------------------------------------------------------------
            j = nsubsplus
            DO i=1,nsubsplus
               IF (min_sub(1,i).EQ.min_phys(1)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                    
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) + len_phys(1)
                  min_sub(2,j) = min_sub(2,i)
                  max_sub(1,j) = max_sub(1,i) + len_phys(1)
                  max_sub(2,j) = max_sub(2,i)
               ENDIF 
               IF (max_sub(1,i).EQ.max_phys(1)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) - len_phys(1)
                  min_sub(2,j) = min_sub(2,i)
                  max_sub(1,j) = max_sub(1,i) - len_phys(1)
                  max_sub(2,j) = max_sub(2,i)
               ENDIF 
            ENDDO
            nsubsplus = j
         ELSE
            !-------------------------------------------------------------------
            !  three dimensions
            !-------------------------------------------------------------------
            j = nsubsplus
            DO i=1,nsubsplus
               IF (min_sub(1,i).EQ.min_phys(1)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) + len_phys(1)
                  min_sub(2,j) = min_sub(2,i)
                  min_sub(3,j) = min_sub(3,i)
                  max_sub(1,j) = max_sub(1,i) + len_phys(1)
                  max_sub(2,j) = max_sub(2,i)
                  max_sub(3,j) = max_sub(3,i)
               ENDIF 
               IF (max_sub(1,i).EQ.max_phys(1)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) - len_phys(1)
                  min_sub(2,j) = min_sub(2,i)
                  min_sub(3,j) = min_sub(3,i)
                  max_sub(1,j) = max_sub(1,i) - len_phys(1)
                  max_sub(2,j) = max_sub(2,i)
                  max_sub(3,j) = max_sub(3,i)
               ENDIF 
            ENDDO
            nsubsplus = j
         ENDIF 
      ENDIF 

      IF (bcdef(3).EQ.ppm_param_bcdef_periodic) THEN
         !----------------------------------------------------------------------
         !  Copy in the second direction
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            !-------------------------------------------------------------------
            !  Two dimensions
            !-------------------------------------------------------------------
            j = nsubsplus
            DO i=1,nsubsplus
               IF (min_sub(2,i).EQ.min_phys(2)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i)
                  min_sub(2,j) = min_sub(2,i) + len_phys(2)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i) + len_phys(2)
               ENDIF 
               IF (max_sub(2,i).EQ.max_phys(2)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) 
                  min_sub(2,j) = min_sub(2,i) - len_phys(2)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i) - len_phys(2)
               ENDIF 
            ENDDO
            nsubsplus = j
         ELSE
            !-------------------------------------------------------------------
            !  three dimensions
            !-------------------------------------------------------------------
            j = nsubsplus 
            DO i=1,nsubsplus
               IF (min_sub(2,i).EQ.min_phys(2)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i)
                  min_sub(2,j) = min_sub(2,i) + len_phys(2)
                  min_sub(3,j) = min_sub(3,i)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i) + len_phys(2)
                  max_sub(3,j) = max_sub(3,i)
               ENDIF 
               IF (max_sub(2,i).EQ.max_phys(2)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) 
                  min_sub(2,j) = min_sub(2,i) - len_phys(2)
                  min_sub(3,j) = min_sub(3,i)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i) - len_phys(2)
                  max_sub(3,j) = max_sub(3,i)
               ENDIF 
            ENDDO
            nsubsplus = j
         ENDIF 
      ENDIF 

      IF (ppm_dim.EQ.3) THEN
         IF (bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
            !-------------------------------------------------------------------
            !  Copy in the third direction
            !-------------------------------------------------------------------
            j = nsubsplus
            DO i=1,nsubsplus
               IF (min_sub(3,i).EQ.min_phys(3)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i)
                  min_sub(2,j) = min_sub(2,i)
                  min_sub(3,j) = min_sub(3,i) + len_phys(3)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i)
                  max_sub(3,j) = max_sub(3,i) + len_phys(3)
               ENDIF 
               IF (max_sub(3,i).EQ.max_phys(3)) THEN
                  j            = j + 1
                  IF (j.GT.isize) THEN
                     isize  = isize * 2
                     ldc(1) = isize
                     CALL ppm_alloc(subid,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     ldc(1) = ppm_dim
                     ldc(2) = isize
                     CALL ppm_alloc(min_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                     CALL ppm_alloc(max_sub,ldc,iopt,info)
                     IF (info.NE.0) GOTO 100
                  ENDIF 
                  subid(j)     = subid(i)
                  min_sub(1,j) = min_sub(1,i) 
                  min_sub(2,j) = min_sub(2,i) 
                  min_sub(3,j) = min_sub(3,i) - len_phys(3)
                  max_sub(1,j) = max_sub(1,i)
                  max_sub(2,j) = max_sub(2,i)
                  max_sub(3,j) = max_sub(3,i) - len_phys(3)
               ENDIF 
            ENDDO
            nsubsplus = j
         ENDIF 
      ENDIF

      !-------------------------------------------------------------------------
      !  Catch allocation errors
      !-------------------------------------------------------------------------
  100 CONTINUE 
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_copy_image_subs',     &
     &        'insufficient memory for sub domains ',__LINE__,info)
          GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_copy_image_subs',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_copy_image_subs_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_copy_image_subs_d
#endif
