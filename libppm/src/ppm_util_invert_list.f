      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_invert_list
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Creates the inverse list (indiced become entries
      !                 and vice versa) for a 1d integer list.
      !
      !  Input        : inlist     (I) original list
      !
      !  Output       : outlist    (I) inverted list
      !                 info       (I) return status. 
      !
      !  Remarks      : Entries in outlist for which no corresponding entry
      !                 in inlist exists are initialized to
      !                 -HUGE(outlist(.)) and can be recognized like this.
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_invert_list.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2005/04/28 12:26:24  pchatela
      !  "Bugfix": Allow the treatment of an empty list -> return an empty list
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_util_invert_list(inlist,outlist,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! non-negative integer list to be inverted
      INTEGER, DIMENSION(:) , POINTER       :: inlist
      ! inverse list
      INTEGER, DIMENSION(:) , POINTER       :: outlist
      ! return status
      INTEGER               , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(1)                 :: ldl,ldu
      INTEGER                               :: iopt
      INTEGER                               :: i,j,inmin,inmax,outmin,outmax
      CHARACTER(LEN=ppm_char)               :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_invert_list',t0,info)
      inmin = LBOUND(inlist,1)
      inmax = UBOUND(inlist,1)
      IF (ppm_debug .GE. 2) THEN
          WRITE(mesg,'(A,2I10)') 'Input list bounds ', inmin, inmax
          CALL ppm_write(ppm_rank,'ppm_util_invert_list',  &
    &                    mesg,info)
      ENDIF
      !-------------------------------------------------------------------------
      !  Determine min and max value in inlist
      !-------------------------------------------------------------------------
      outmax = -HUGE(outmax)
      outmin = HUGE(outmin)
      DO i=inmin,inmax
          IF (inlist(i) .GT. outmax) outmax = inlist(i)
          IF (inlist(i) .LT. outmin) outmin = inlist(i)
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Case of an empty input list
      !-------------------------------------------------------------------------
      IF (outmax .LT. outmin) THEN
          ! this is equivalent to testing (inmax .LT. inmin)
          ! We allocate an empty inverse list
          iopt = ppm_param_alloc_fit
          ldu(1) = 0
          CALL ppm_alloc(outlist,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_invert_list',        &
    &              'inverted list OUTLIST',__LINE__,info)
              GOTO 9999
          ENDIF
      !-------------------------------------------------------------------------
      !  General case: non-empty list
      !-------------------------------------------------------------------------
      ELSE
      !-------------------------------------------------------------------------
      !  Allocate outlist. No need to preserve its contents
      !-------------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldl(1) = outmin
          ldu(1) = outmax
          IF (ppm_debug .GE. 2) THEN
              WRITE(mesg,'(A,2I10)') 'Output list bounds ', outmin, outmax
              CALL ppm_write(ppm_rank,'ppm_util_invert_list',  &
    &              mesg,info)
          ENDIF
          CALL ppm_alloc(outlist,ldl,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_invert_list',  &
    &                        'inverted list OUTLIST',__LINE__,info)
              GOTO 9999
          ENDIF

      !-------------------------------------------------------------------------
      !  Initialize outlist
      !-------------------------------------------------------------------------
          DO i=outmin,outmax
              outlist(i) = -HUGE(outlist(i))
          ENDDO

      !-------------------------------------------------------------------------
      !  Build inverse list
      !-------------------------------------------------------------------------
          DO i=inmin,inmax
              outlist(inlist(i)) = i
          ENDDO
      END IF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_invert_list',t0,info)
      RETURN
      END SUBROUTINE ppm_util_invert_list
