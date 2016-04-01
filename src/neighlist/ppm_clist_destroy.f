      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_clist_destroy
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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

      SUBROUTINE ppm_clist_destroy(clist,info)
      !!! Properly deallocates the cell list it is passed.
      !!!
      !!! [NOTE]
      !!! At least using pgf90, this routine is actually not necessary as
      !!! `DEALLOCATE(clist)` would be sufficient (no memory leak would
      !!! occur according to Valgrind). But since this might be a
      !!! compiler-dependent feature we do it the orthodox way for the sake
      !!! of portability.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_t_clist), DIMENSION(:), POINTER       :: clist
      !!! Cell list which is to be deallocated.
      INTEGER,                         INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double) :: t0
      ! counters
      INTEGER               :: i
      ! for allocate
      INTEGER, DIMENSION(2) :: lda
      INTEGER               :: iopt

      CHARACTER(LEN=ppm_char) :: caller='ppm_clist_destroy'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Free cell list memory. At least using pgf90, this is actually not
      !  necessary as DEALLOCATE(clist) would be sufficient (no memory leak
      !  would occur according to Valgrind). But since this might be a
      !  compiler-dependent feature we do it the orthodox way for the sake
      !  of portability.
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ASSOCIATED(clist)) THEN
         DO i=1,SIZE(clist,1)
            CALL ppm_alloc(clist(i)%nm,lda,iopt,info)
            or_fail_dealloc('cell list CLIST%LHBX')

            CALL ppm_alloc(clist(i)%lhbx,lda,iopt,info)
            or_fail_dealloc('cell list CLIST%NM')

            CALL ppm_alloc(clist(i)%lpdx,lda,iopt,info)
            or_fail_dealloc('cell list CLIST%LPDX')
         ENDDO
         DEALLOCATE(clist, STAT=info)
         or_fail_dealloc('cell list CLIST')
         NULLIFY(clist)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
             fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
          ENDIF
          IF (.NOT. ASSOCIATED(clist)) THEN
             fail('clist is not associated!',exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_clist_destroy
