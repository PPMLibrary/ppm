      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_copy
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

      SUBROUTINE ppm_topo_copy(intopo,outtopo,info)
      !!! This routine copies the member structures from one topology
      !!! object to another one. The previous contents of outtopo will
      !!! be destroyed. This also allocates outtopo.
      !!!
      !!! [WARNING]
      !!! This routine is not tested

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_typedef
      USE ppm_module_topo_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), INTENT(IN   ) :: intopo
      !!! Input topology
      TYPE(ppm_t_topo), INTENT(  OUT) :: outtopo
      !!! Output topology. Will contain the same data as intopo upon return.
      INTEGER,          INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)  :: t0

      INTEGER                :: nsubs,nsublist,maxneigh,prec,i,j,iopt
      INTEGER, DIMENSION(1)  :: ldc

      CHARACTER(LEN=ppm_char) :: caller='ppm_topo_copy'

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
      !  Determine size of intopo
      !-------------------------------------------------------------------------
      nsubs    = intopo%nsubs
      nsublist = intopo%nsublist
      maxneigh = SIZE(intopo%ineighsubs,1)
      prec     = intopo%prec

      !-------------------------------------------------------------------------
      !  Allocate result topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_alloc(outtopo%ID,nsubs,nsublist,maxneigh,prec,info)
      IF (info .NE. ppm_param_success) THEN
          fail('result topology OUTTOPO',ppm_err_alloc,ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the scalar members
      !-------------------------------------------------------------------------
      outtopo%ID          = intopo%ID
      outtopo%prec        = intopo%prec
      outtopo%nsubs       = intopo%nsubs
      outtopo%nsublist    = intopo%nsublist
      outtopo%nneighproc  = intopo%nneighproc
      outtopo%isoptimized = intopo%isoptimized
      outtopo%ncommseq    = intopo%ncommseq

      !-------------------------------------------------------------------------
      !  Copy the boundary conditions
      !-------------------------------------------------------------------------
      DO i=1,2*ppm_dim
          outtopo%bcdef(i) = intopo%bcdef(i)
      ENDDO

      SELECT CASE (prec)
      CASE (ppm_kind_single)
          SELECT CASE (ppm_dim)
          CASE (2)
             !------------------------------------------------------------------
             !  Copy the domain
             !------------------------------------------------------------------
             outtopo%min_physs(1) = intopo%min_physs(1)
             outtopo%min_physs(2) = intopo%min_physs(2)
             outtopo%max_physs(1) = intopo%max_physs(1)
             outtopo%max_physs(2) = intopo%max_physs(2)
             !------------------------------------------------------------------
             !  Copy the subdomains
             !------------------------------------------------------------------
             DO i=1,nsubs
                outtopo%min_subs(1,i) = intopo%min_subs(1,i)
                outtopo%min_subs(2,i) = intopo%min_subs(2,i)
                outtopo%max_subs(1,i) = intopo%max_subs(1,i)
                outtopo%max_subs(2,i) = intopo%max_subs(2,i)
             ENDDO

          CASE DEFAULT
             !------------------------------------------------------------------
             !  Copy the domain
             !------------------------------------------------------------------
             outtopo%min_physs(1) = intopo%min_physs(1)
             outtopo%min_physs(2) = intopo%min_physs(2)
             outtopo%min_physs(3) = intopo%min_physs(3)
             outtopo%max_physs(1) = intopo%max_physs(1)
             outtopo%max_physs(2) = intopo%max_physs(2)
             outtopo%max_physs(3) = intopo%max_physs(3)
             !------------------------------------------------------------------
             !  Copy the subdomains
             !------------------------------------------------------------------
             DO i=1,nsubs
                outtopo%min_subs(1,i) = intopo%min_subs(1,i)
                outtopo%min_subs(2,i) = intopo%min_subs(2,i)
                outtopo%min_subs(3,i) = intopo%min_subs(3,i)
                outtopo%max_subs(1,i) = intopo%max_subs(1,i)
                outtopo%max_subs(2,i) = intopo%max_subs(2,i)
                outtopo%max_subs(3,i) = intopo%max_subs(3,i)
             ENDDO

          END SELECT
          !---------------------------------------------------------------------
          !  Copy the sub costs
          !---------------------------------------------------------------------
          DO i=1,nsubs
             outtopo%sub_costs(i) = intopo%sub_costs(i)
          ENDDO

      CASE DEFAULT
          SELECT CASE (ppm_dim)
          CASE (2)
             !------------------------------------------------------------------
             !  Copy the domain
             !------------------------------------------------------------------
             outtopo%min_physd(1) = intopo%min_physd(1)
             outtopo%min_physd(2) = intopo%min_physd(2)
             outtopo%max_physd(1) = intopo%max_physd(1)
             outtopo%max_physd(2) = intopo%max_physd(2)
             !------------------------------------------------------------------
             !  Copy the subdomains
             !------------------------------------------------------------------
             DO i=1,nsubs
                outtopo%min_subd(1,i) = intopo%min_subd(1,i)
                outtopo%min_subd(2,i) = intopo%min_subd(2,i)
                outtopo%max_subd(1,i) = intopo%max_subd(1,i)
                outtopo%max_subd(2,i) = intopo%max_subd(2,i)
             ENDDO

          CASE DEFAULT
             !------------------------------------------------------------------
             !  Copy the domain
             !------------------------------------------------------------------
             outtopo%min_physd(1) = intopo%min_physd(1)
             outtopo%min_physd(2) = intopo%min_physd(2)
             outtopo%min_physd(3) = intopo%min_physd(3)
             outtopo%max_physd(1) = intopo%max_physd(1)
             outtopo%max_physd(2) = intopo%max_physd(2)
             outtopo%max_physd(3) = intopo%max_physd(3)
             !------------------------------------------------------------------
             !  Copy the subdomains
             !------------------------------------------------------------------
             DO i=1,nsubs
                outtopo%min_subd(1,i) = intopo%min_subd(1,i)
                outtopo%min_subd(2,i) = intopo%min_subd(2,i)
                outtopo%min_subd(3,i) = intopo%min_subd(3,i)
                outtopo%max_subd(1,i) = intopo%max_subd(1,i)
                outtopo%max_subd(2,i) = intopo%max_subd(2,i)
                outtopo%max_subd(3,i) = intopo%max_subd(3,i)
             ENDDO

          END SELECT
          !---------------------------------------------------------------------
          !  Copy the sub costs
          !---------------------------------------------------------------------
          DO i=1,nsubs
              outtopo%sub_costd(i) = intopo%sub_costd(i)
          ENDDO

      END SELECT

      !-------------------------------------------------------------------------
      !  Copy the sub to proc assignment
      !-------------------------------------------------------------------------
      DO i=1,nsubs
          outtopo%sub2proc(i) = intopo%sub2proc(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy the list of local subs
      !-------------------------------------------------------------------------
      DO i=1,nsublist
          outtopo%isublist(i) = intopo%isublist(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy the subs boundary conditions
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO i=1,nsublist
              outtopo%subs_bc(1,i) = intopo%subs_bc(1,i)
              outtopo%subs_bc(2,i) = intopo%subs_bc(2,i)
              outtopo%subs_bc(3,i) = intopo%subs_bc(3,i)
              outtopo%subs_bc(4,i) = intopo%subs_bc(4,i)
          ENDDO
      ELSE
          DO i=1,nsublist
              outtopo%subs_bc(1,i) = intopo%subs_bc(1,i)
              outtopo%subs_bc(2,i) = intopo%subs_bc(2,i)
              outtopo%subs_bc(3,i) = intopo%subs_bc(3,i)
              outtopo%subs_bc(4,i) = intopo%subs_bc(4,i)
              outtopo%subs_bc(5,i) = intopo%subs_bc(5,i)
              outtopo%subs_bc(6,i) = intopo%subs_bc(6,i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the sub neighbor numbers
      !-------------------------------------------------------------------------
      DO i=1,nsublist
          outtopo%nneighsubs(i) = intopo%nneighsubs(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy the sub neighbor lists
      !-------------------------------------------------------------------------
      DO i=1,nsublist
          DO j=1,maxneigh
              outtopo%ineighsubs(j,i) = intopo%ineighsubs(j,i)
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy the list of neighboring processors. This is not allocated by
      !  topo_alloc!
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = intopo%nneighproc
      CALL ppm_alloc(outtopo%ineighproc,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
         fail('list of neighboring proc OUTTOPO%INEIGHPROC',ppm_err_alloc,ppm_error=ppm_error_fatal)
      ENDIF
      DO i=1,intopo%nneighproc
         outtopo%ineighproc(i) = intopo%ineighproc(i)
      ENDDO

      IF (intopo%isoptimized) THEN
          iopt   = ppm_param_alloc_fit
          ldc(1) = intopo%ncommseq
          CALL ppm_alloc(outtopo%icommseq,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
             fail('communication sequence OUTTOPO%ICOMMSEQ',ppm_err_alloc,ppm_error=ppm_error_fatal)
          ENDIF
          DO i=1,intopo%ncommseq
             outtopo%icommseq(i) = intopo%icommseq(i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Set outtopo to defined
      !-------------------------------------------------------------------------
      outtopo%isdefined = .TRUE.

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
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_topo_copy
