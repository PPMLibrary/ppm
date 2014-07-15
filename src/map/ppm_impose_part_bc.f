      !-------------------------------------------------------------------------
      !  Subroutine   :           ppm_impose_part_bc
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_impose_part_bc_s(topoid,xp,Npart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_impose_part_bc_d(topoid,xp,Npart,info)
#endif
      !!! This routine imposes periodic boundary conditions for
      !!! the particles. That is it moves/wraps any particles
      !!! outside the computational box into the box. The routine
      !!! assumes that the particles are located within the size
      !!! of one computational box.
      !!!
      !!! [NOTE]
      !!! This routine is a user callable routine.
      !!!
      !!! [CAUTION]
      !!! (This concerns only PPM developers)                                  +
      !!! The routine is sensitive the round off errors! So be careful not to
      !!! change the order of the IF statements (see code comments) and make
      !!! sure the compiler is not doing something wrong as well! Originally
      !!! when particles would be close/at the boundary, the routine would
      !!! fail to map all particles, but succed if called a 2nd time. The
      !!! changed in the order of the IF statements removed this double calling
      !!! sequence.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_data
      USE ppm_module_check_id
      USE ppm_module_topo_typedef

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      !!! Particle coordinates
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      REAL(MK), DIMENSION(ppm_dim) :: len_phys,max_phys,min_phys
      REAL(MK)                     :: t0

      INTEGER :: i

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller='ppm_impose_part_bc'

      LOGICAL :: valid

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  compute the size of the computational domain
      !-------------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
      len_phys(1:ppm_dim)=topo%max_physd(1:ppm_dim)-topo%min_physd(1:ppm_dim)
      max_phys(1:ppm_dim)=topo%max_physd(1:ppm_dim)
      min_phys(1:ppm_dim)=topo%min_physd(1:ppm_dim)
#else
      len_phys(1:ppm_dim)=topo%max_physs(1:ppm_dim)-topo%min_physs(1:ppm_dim)
      max_phys(1:ppm_dim)=topo%max_physs(1:ppm_dim)
      min_phys(1:ppm_dim)=topo%min_physs(1:ppm_dim)
#endif

      !-------------------------------------------------------------------------
      !  wrap the particles on this topology - do that by allowing particles to
      !  stay on the left boundary but always wrap particle located on the right
      !  boundary. This GE and LT MUST be consistent with the logic of the cell
      !  lists to allow the FLOOR in these routines to NOT to return Nm(:)+1,
      !  which would be the case if the particle was allowed to stay on the
      !  right boundary.
      !-------------------------------------------------------------------------
      IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN
         DO i=1,Npart
            !-------------------------------------------------------------------
            !  The order of the if statements matter !
            !  and you need both if statements separately (not ELSEIFs) !
            !  all this due to round off errors
            !-------------------------------------------------------------------
            IF (xp(1,i).LT.min_phys(1)) THEN
               xp(1,i) = xp(1,i) + len_phys(1)
            ENDIF
            IF (xp(1,i).GE.max_phys(1)) THEN
               xp(1,i) = xp(1,i) - len_phys(1)
            ENDIF
         ENDDO
      ENDIF
      IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN
         DO i=1,Npart
            IF (xp(2,i).LT.min_phys(2)) THEN
               xp(2,i) = xp(2,i) + len_phys(2)
            ENDIF
            IF (xp(2,i).GE.max_phys(2)) THEN
               xp(2,i) = xp(2,i) - len_phys(2)
            ENDIF
         ENDDO
      ENDIF
      IF (ppm_dim.EQ.3) THEN
         IF (topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
            DO i=1,Npart
               IF (xp(3,i).LT.min_phys(3)) THEN
                  xp(3,i) = xp(3,i) + len_phys(3)
               ENDIF
               IF (xp(3,i).GE.max_phys(3)) THEN
                  xp(3,i) = xp(3,i) - len_phys(3)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        CALL ppm_check_topoid(topoid,valid,info)
        IF (.NOT. valid) THEN
           fail('topoid is invalid!',exit_point=8888)
        ENDIF
        IF (Npart .LT. 0) THEN
           fail('Npart must be >= 0',exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_impose_part_bc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_impose_part_bc_d
#endif
