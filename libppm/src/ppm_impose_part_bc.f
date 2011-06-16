      !-------------------------------------------------------------------------
      !  Subroutine   :           ppm_impose_part_bc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine imposes periodic boundary conditions for 
      !                 the particles. That is it moves/wraps any particles 
      !                 outside the computational box into the box. The routine
      !                 assumes that the particles are located within the size
      !                 of one computational box.
      !  
      !  Input        : Npart        (I) the number of particles
      !                 topo_id      (I) user topology (not ppm internal!).
      !
      !  Input/output : xp(:.:)      (F) particle position
      !
      !  Output       : info         (I) return status: zero on success
      !
      !  Routines     : ppm_alloc
      !                 ppm_error
      !                 substart
      !                 substop
      !
      !  Remarks      : This routine is a user callable routine
      !                 The routine is sensitive the round off errors !
      !                 So be careful not to change the order of the IF
      !                 statements (see below) and make sure the compiler
      !                 is not doing something wrong as well! Originally
      !                 when particles would be close/at the boundary, the
      !                 routine would fail to map all particles, but succed
      !                 if called a 2nd time. The changed in the order of the
      !                 IF statements removed this double calling sequence.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_impose_part_bc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/09/04 18:34:47  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.4  2004/09/28 21:42:16  walther
      !  Bug fix: due to round off, the order of the IF statements matter.
      !
      !  Revision 1.3  2004/08/31 13:29:57  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.2  2004/08/03 11:38:32  ivos
      !  bugfix: argument check was on ppm_topoid instead of topo_id.
      !  Fixed wrong output message. Fixed wrong indices in ppm_bcdef.
      !  Routine tested in 2d and 3d.
      !
      !  Revision 1.1  2004/08/03 08:07:10  walther
      !  Initial version (untested).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_impose_part_bc_s(xp,Npart,topo_id,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_impose_part_bc_d(xp,Npart,topo_id,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_data
      USE ppm_module_check_topoid

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
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
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      INTEGER                 , INTENT(IN   ) :: topo_id ! user topoid
      INTEGER                 , INTENT(IN   ) :: Npart 
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                     :: t0
      INTEGER                      :: topoid ! ppm internal topoid
      INTEGER                      :: i
      CHARACTER(LEN=ppm_char)      :: mesg
      LOGICAL                      :: valid
      REAL(MK), DIMENSION(ppm_dim) :: len_phys,max_phys,min_phys
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_impose_part_bc',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
         IF (.NOT. valid) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_impose_part_bc',  &
     &           'topo_id is invalid!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (Npart .LT. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_impose_part_bc',  &
     &           'Npart must be >= 0',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  find the corresponding internal ppm topoid (topoid)
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)

      !-------------------------------------------------------------------------
      !  compute the size of the computational domain
      !-------------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
      len_phys(1:ppm_dim) = ppm_max_physd(1:ppm_dim,topoid) -    &
     &                      ppm_min_physd(1:ppm_dim,topoid)
      max_phys(1:ppm_dim) = ppm_max_physd(1:ppm_dim,topoid)
      min_phys(1:ppm_dim) = ppm_min_physd(1:ppm_dim,topoid)
#else
      len_phys(1:ppm_dim) = ppm_max_physs(1:ppm_dim,topoid) -    &
     &                      ppm_min_physs(1:ppm_dim,topoid)
      max_phys(1:ppm_dim) = ppm_max_physs(1:ppm_dim,topoid)
      min_phys(1:ppm_dim) = ppm_min_physs(1:ppm_dim,topoid)
#endif

      !-------------------------------------------------------------------------
      !  wrap the particles on this topology - do that by allowing particles to
      !  stay on the left boundary but always wrap particle located on the right
      !  boundary. This GE and LT MUST be consistent with the logic of the cell
      !  lists to allow the FLOOR in these routines to NOT to return Nm(:)+1,
      !  which would be the case if the particle was allowed to stay on the 
      !  right boundary.
      !-------------------------------------------------------------------------
      IF (ppm_bcdef(1,topoid).EQ.ppm_param_bcdef_periodic) THEN
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
      IF (ppm_bcdef(3,topoid).EQ.ppm_param_bcdef_periodic) THEN
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
         IF (ppm_bcdef(5,topoid).EQ.ppm_param_bcdef_periodic) THEN
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
      CALL substop('ppm_impose_part_bc',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_impose_part_bc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_impose_part_bc_d
#endif
