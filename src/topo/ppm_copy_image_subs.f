      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_copy_image_subs
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE copy_imgsubs_s(min_phys,max_phys,bcdef, &
     &   min_sub,max_sub,nsubs,subid,nsubsplus,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE copy_imgsubs_d(min_phys,max_phys,bcdef, &
     &   min_sub,max_sub,nsubs,subid,nsubsplus,info)
#endif
      !!! This routine copies the sub domains to ghost sub domains
      !!! in the case of periodic boundary condition and stores
      !!! the original ID of the sub domain - this to allow the
      !!! find neighbours to search (without wrapping space) the
      !!! list of sub domains and still retrieve the real ID of sub domain.
      !!!
      !!! [NOTE]
      !!! This is a happy going routine: we are comparing the equality of
      !!! floating point numbers *but* it *works* since the value of subs
      !!! are *assigned* from the phys min/max and *not* computed (byte copy).
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys
      !!! Min. extent of the physical domain
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: max_phys
      !!! Max. extent of the physical domain
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: bcdef
      !!! Boundary condition definition
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Min. extent of the sub domain including after the call to this
      !!! routine the ghost sub domains
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Max. extent of the sub domain including after the call to this
      !!! routine the ghost sub domains
      INTEGER , DIMENSION(:)  , POINTER       :: subid
      !!! Real id of the sub domain
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Total number of sub domains
      INTEGER                 , INTENT(  OUT) :: nsubsplus
      !!! Total number of sub domains plus the ghost sub domains
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      REAL(MK)                     :: t0
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER                      :: i,j,k
      INTEGER                      :: istat,iopt,isize
      CHARACTER(ppm_char)          :: mesg
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
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
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
      CONTAINS
       SUBROUTINE check
            IF (.NOT. ASSOCIATED(min_sub)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'min_sub should be associated!',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (.NOT. ASSOCIATED(max_sub)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'max_sub should be associated!',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (.NOT. ASSOCIATED(subid)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'subid should be associated!',__LINE__,info)
                GOTO 8888
            ENDIF
        !----------------------------------------------------------------------
        !  The physical boundary is sane ?
        !----------------------------------------------------------------------
            IF (max_phys(1) .LE. min_phys(1)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'max_phys(1) must be > min_phys(1)',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (max_phys(2) .LE. min_phys(2)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'max_phys(2) must be > min_phys(2)',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (ppm_dim .GT. 2) THEN
                IF (max_phys(3) .LE. min_phys(3)) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',&
     &                     'max_phys(3) must be > min_phys(3)',__LINE__,info)
                GOTO 8888
                ENDIF
            ENDIF

        !----------------------------------------------------------------------
        !  The periodic BC are sane ?
        !----------------------------------------------------------------------
            IF (bcdef(1).EQ.ppm_param_bcdef_periodic.AND. &
     &          bcdef(2).NE.ppm_param_bcdef_periodic) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'periodic on only one x face !',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (bcdef(3).EQ.ppm_param_bcdef_periodic.AND. &
     &          bcdef(4).NE.ppm_param_bcdef_periodic) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                     'periodic on only one y face !',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (ppm_dim.GT.2) THEN
                IF (bcdef(5).EQ.ppm_param_bcdef_periodic.AND. &
     &              bcdef(6).NE.ppm_param_bcdef_periodic) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_copy_image_subs',  &
     &                        'periodic on only one z face !',__LINE__,info)
                GOTO 8888
                ENDIF
            ENDIF
 8888       CONTINUE
        END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE copy_imgsubs_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE copy_imgsubs_d
#endif
