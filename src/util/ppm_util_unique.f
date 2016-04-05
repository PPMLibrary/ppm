      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_unique
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 CSE Lab (ETH Zurich), MOSAIC Group (MPI-CBG Dresden),
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_unique_s(inlist,outlist,info,inlistSize,outlistsize)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_unique_d(inlist,outlist,info,inlistSize,outlistsize)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_util_unique_i(inlist,outlist,info,inlistSize,outlistsize)
#elif __KIND == __LONGINT
      SUBROUTINE ppm_util_unique_li(inlist,outlist,info,inlistSize,outlistSize)
#endif
      !!! Taking the input List of values, create the sorted unique list
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_qsort, ONLY : ppm_util_qsort
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION || __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: inlist
      !!! input List
      REAL(MK), DIMENSION(:), POINTER       :: outlist
      !!! sorted unique output List
#elif __KIND == __INTEGER
      INTEGER,  DIMENSION(:), INTENT(IN   ) :: inlist
      !!! input List
      INTEGER,  DIMENSION(:), POINTER       :: outlist
      !!! sorted unique output List
#elif __KIND == __LONGINT
      INTEGER(ppm_kind_int64), DIMENSION(:), INTENT(IN   ) :: inlist
      !!! input List
      INTEGER(ppm_kind_int64), DIMENSION(:), POINTER       :: outlist
      !!! sorted unique output List
#endif
      INTEGER,                INTENT(OUT)   :: info
      !!! Return status, 0 on success
      INTEGER, OPTIONAL,      INTENT(IN   ) :: inlistSize
      !!! Size of the elements in the input array to make unique
      INTEGER, OPTIONAL,      INTENT(  OUT) :: outlistSize
      !!! size of the output unique elements

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
#if   __KIND == __SINGLE_PRECISION || __KIND == __DOUBLE_PRECISION
      REAL(MK)              :: ppm_myeps
#endif

      INTEGER                         :: inSize
      INTEGER                         :: outSize
      INTEGER,  DIMENSION(:), POINTER :: indxlist
      INTEGER, DIMENSION(1)           :: ldu
      INTEGER                         :: iopt
      INTEGER                         :: i
      INTEGER                         :: j

      CHARACTER(LEN=ppm_char) :: caller='ppm_util_unique'

      LOGICAL, DIMENSION(SIZE(inlist)) :: MASK
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialization
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      IF (PRESENT(inlistSize)) THEN
         inSize=inlistSize
      ELSE
         inSize=SIZE(inlist,DIM=1)
      ENDIF
      IF (inSize.EQ.0) GOTO 9999

#if   __KIND == __SINGLE_PRECISION
      ppm_myeps=ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      ppm_myeps=ppm_myepsd
#endif

      NULLIFY(indxlist)
      CALL ppm_util_qsort(inlist,indxlist,info,inSize)
      or_fail("ppm_util_qsort")

      outSize=1
      DO i=2,inSize
#if   __KIND == __SINGLE_PRECISION || __KIND == __DOUBLE_PRECISION
         ! If the difference between two elements is less than Floating point tolerance
         ! then I would consider them as equal
         IF ((inlist(indxlist(i))-inlist(indxlist(i-1))).LE.ppm_myeps) THEN
#elif __KIND == __INTEGER || __KIND == __LONGINT
         IF (inlist(indxlist(i)).EQ.inlist(indxlist(i-1))) THEN
#endif
            MASK(indxlist(i))=.FALSE.
         ELSE
            MASK(indxlist(i))=.TRUE.
            outSize=outSize+1
         ENDIF
      ENDDO

      iopt=ppm_param_alloc_fit
      ldu=outSize
      CALL ppm_alloc(outlist,ldu,iopt,info)
      or_fail_alloc("outlist")

      outlist(1)=inlist(indxlist(1))
      j=1
      DO i=2,inSize
         IF (MASK(i)) THEN
            j=j+1
            outlist(j)=inlist(indxlist(i))
         ENDIF
      ENDDO

      ! Free memory
      iopt=ppm_param_dealloc
      CALL ppm_alloc(indxlist,ldu,iopt,info)
      or_fail_dealloc("indxlist")

      IF (PRESENT(outlistSize)) outlistSize=outSize

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_unique_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_unique_d
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_util_unique_i
#elif __KIND == __LONGINT
      END SUBROUTINE ppm_util_unique_li
#endif