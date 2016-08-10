      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_1d
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


      SUBROUTINE ppm_alloc_topo(topo,iopt,info)
      !!! (Re)allocates the memory of one-dimensional arrays
      !!! (pointers) based on the number of elements.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_topo_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER       :: topo
      INTEGER,          INTENT(IN   ) :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_dealloc
      INTEGER,          INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------

      LOGICAL :: lalloc,ldealloc

      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=*), PARAMETER :: caller='ppm_alloc_topo'

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
      ! maybe add some sanity checks later

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lalloc   = .FALSE.
      ldealloc = .FALSE.

      SELECT CASE (iopt)
      CASE (ppm_param_alloc_fit)
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(topo)) ldealloc = .TRUE.
         lalloc   = .TRUE.

      CASE (ppm_param_dealloc)
         ldealloc = .TRUE.

      CASE DEFAULT
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         fail('unknown iopt')

      END SELECT

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate old data first
      !-------------------------------------------------------------------------
      IF (ldealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(topo)) THEN
            ! first deallocate all content of topo
            IF (ASSOCIATED(topo%min_physs))   DEALLOCATE(topo%min_physs,  STAT=info)
            IF (ASSOCIATED(topo%max_physs))   DEALLOCATE(topo%max_physs,  STAT=info)
            IF (ASSOCIATED(topo%min_physd))   DEALLOCATE(topo%min_physd,  STAT=info)
            IF (ASSOCIATED(topo%max_physd))   DEALLOCATE(topo%max_physd,  STAT=info)
            IF (ASSOCIATED(topo%bcdef))       DEALLOCATE(topo%bcdef,      STAT=info)
            IF (ASSOCIATED(topo%min_subs))    DEALLOCATE(topo%min_subs,   STAT=info)
            IF (ASSOCIATED(topo%max_subs))    DEALLOCATE(topo%max_subs,   STAT=info)
            IF (ASSOCIATED(topo%min_subd))    DEALLOCATE(topo%min_subd,   STAT=info)
            IF (ASSOCIATED(topo%max_subd))    DEALLOCATE(topo%max_subd,   STAT=info)
            IF (ASSOCIATED(topo%sub_costs))   DEALLOCATE(topo%sub_costs,  STAT=info)
            IF (ASSOCIATED(topo%sub_costd))   DEALLOCATE(topo%sub_costd,  STAT=info)
            IF (ASSOCIATED(topo%sub2proc))    DEALLOCATE(topo%sub2proc,   STAT=info)
            IF (ASSOCIATED(topo%isublist))    DEALLOCATE(topo%isublist,   STAT=info)
            IF (ASSOCIATED(topo%subs_bc))     DEALLOCATE(topo%subs_bc,    STAT=info)
            IF (ASSOCIATED(topo%ineighsubs))  DEALLOCATE(topo%ineighsubs, STAT=info)
            IF (ASSOCIATED(topo%nneighsubs))  DEALLOCATE(topo%nneighsubs, STAT=info)
            IF (ASSOCIATED(topo%ineighproc))  DEALLOCATE(topo%ineighproc, STAT=info)
            IF (ASSOCIATED(topo%ineighcolor)) DEALLOCATE(topo%ineighcolor,STAT=info)
            IF (ASSOCIATED(topo%icommseq))    DEALLOCATE(topo%icommseq,   STAT=info)
            or_fail_dealloc('Deallocating topo',exit_point=no)

            NULLIFY(topo%min_physs)
            NULLIFY(topo%min_physd)
            NULLIFY(topo%max_physs)
            NULLIFY(topo%max_physd)
            NULLIFY(topo%bcdef)
            NULLIFY(topo%min_subs)
            NULLIFY(topo%max_subs)
            NULLIFY(topo%min_subd)
            NULLIFY(topo%max_subd)
            NULLIFY(topo%sub_costs)
            NULLIFY(topo%sub_costd)
            NULLIFY(topo%sub2proc)
            NULLIFY(topo%isublist)
            NULLIFY(topo%subs_bc)
            NULLIFY(topo%ineighsubs)
            NULLIFY(topo%nneighsubs)
            NULLIFY(topo%ineighproc)
            NULLIFY(topo%ineighcolor)
            NULLIFY(topo%icommseq)

            topo%ID         = ppm_param_undefined
            topo%prec       = ppm_param_undefined
            topo%nsubs      = ppm_param_undefined
            topo%nsublist   = ppm_param_undefined
            topo%nneighproc = ppm_param_undefined
            topo%ncommseq   = ppm_param_undefined

            topo%ghostsizes = 0.0_ppm_kind_single
            topo%ghostsized = 0.0_ppm_kind_double

            DEALLOCATE(topo,STAT=info)
            or_fail_dealloc('Deallocating topo',exit_point=no)

            NULLIFY(topo)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(topo,STAT=info)
         or_fail_alloc('Allocating topo',ppm_error=ppm_error_fatal)

         NULLIFY(topo%min_physs)
         NULLIFY(topo%min_physd)
         NULLIFY(topo%max_physs)
         NULLIFY(topo%max_physd)
         NULLIFY(topo%bcdef)
         NULLIFY(topo%min_subs)
         NULLIFY(topo%max_subs)
         NULLIFY(topo%min_subd)
         NULLIFY(topo%max_subd)
         NULLIFY(topo%sub_costs)
         NULLIFY(topo%sub_costd)
         NULLIFY(topo%sub2proc)
         NULLIFY(topo%isublist)
         NULLIFY(topo%subs_bc)
         NULLIFY(topo%ineighsubs)
         NULLIFY(topo%nneighsubs)
         NULLIFY(topo%ineighproc)
         NULLIFY(topo%ineighcolor)
         !-------------------------------------------------------------------------
         !  Nullify the communication sequence. It will be allocated when first used.
         !-------------------------------------------------------------------------
         NULLIFY(topo%icommseq)

         topo%ID = 0
         topo%prec = ppm_param_undefined
         topo%nsubs = 0
         topo%nsublist = 0
         topo%nneighproc = 0
         topo%ncommseq = 0
         topo%ghostsizes = 0.0_ppm_kind_single
         topo%ghostsized = 0.0_ppm_kind_double
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE ppm_alloc_topo
