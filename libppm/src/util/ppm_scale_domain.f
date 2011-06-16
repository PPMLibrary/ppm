      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_scale_domain
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_scale_domain_s(topoid,scale,origo,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_scale_domain_d(topoid,scale,origo,info)
#endif
      !!! This routine scales the computational domain ie. it scales the
      !!! `topo%ppm_min_phys` and `topo%ppm_max_phys` and the subs associated
      !!! with the topology, `topo` being the topology at `topoid`. The scaling
      !!! is given by: `x_new = (x_old - x_origo)*x_scale + x_origo`
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_check_id
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
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: scale
      !!! Scale factors
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: origo
      !!! Origin of the scaling
      INTEGER               , INTENT(IN   ) :: topoid
      !!! User topology:
      !!!
      !!! * topoid points to a valid topology: scale this topology
      !!! * topoid == ppm_param_topo_undefined: scale all topologies
      INTEGER               , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
      INTEGER                           :: i,j,k,dim
      CHARACTER(LEN=ppm_char)           :: mesg
      LOGICAL                           :: valid
      TYPE(ppm_t_topo)        , POINTER :: topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_scale_domain',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  let us tell the compiler what we want to do ... ie. we should show it
      !  that the first dimension is 1) short and either 2 or 3.
      !-------------------------------------------------------------------------

      IF (ppm_dim.EQ.2) THEN
         dim = 2
      ELSE
         dim = 3
      ENDIF

      !-------------------------------------------------------------------------
      !  based on the topoid we have to scale all or one topology
      !-------------------------------------------------------------------------
      ! TODO: HERE WE MERGED CODE FROM THE DEM SOURCES - IT NEEDS TESTING
      !-------------------------------------------------------------------------
      IF (topoid.EQ.ppm_param_topo_undefined) THEN
         !----------------------------------------------------------------------
         !  scale all topologies
         !----------------------------------------------------------------------

         DO k=1,SIZE(ppm_topo) ! loop over all valid topologies
            CALL ppm_check_topoid(k,valid,info)
            IF (.NOT. valid) THEN
                CONTINUE
            ENDIF
            topo => ppm_topo(k)%t

            DO j=1,topo%nsubs
                DO i=1,dim
#if __KIND == __DOUBLE_PRECISION
                    topo%min_subd(i,j) = &
     &                (topo%min_subd(i,j) - origo(i))*scale(i) + origo(i)
                    topo%max_subd(i,j) = &
     &                (topo%max_subd(i,j) - origo(i))*scale(i) + origo(i)
#else
                    topo%min_subs(i,j) = &
     &                (topo%min_subs(i,j) - origo(i))*scale(i) + origo(i)
                    topo%max_subs(i,j) = &
     &                (topo%max_subs(i,j) - origo(i))*scale(i) + origo(i)
#endif
                ENDDO
            ENDDO ! end loop of subs in k-th topo

            !-------------------------------------------------------------------
            !  scale the physical domain for/in each topology
            !-------------------------------------------------------------------

	        DO i=1,dim
#if __KIND == __DOUBLE_PRECISION
				topo%min_physd(i) = &
     &            (topo%min_physd(i) - origo(i))*scale(i) + origo(i)
                topo%max_physd(i) = &
     &            (topo%max_physd(i) - origo(i))*scale(i) + origo(i)
#else
               topo%min_physs(i) = &
     &           (topo%min_physs(i) - origo(i))*scale(i) + origo(i)
               topo%max_physs(i) = &
     &           (topo%max_physs(i) - origo(i))*scale(i) + origo(i)
#endif
            ENDDO
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  else, scale toplogy with given id
         !----------------------------------------------------------------------
         CALL ppm_check_topoid(topoid,valid,info)
         IF (.NOT. valid) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_scale_domain',  &
     &                     'Given topoid is not valid',__LINE__,info)
            GOTO 9999
         ENDIF
         topo => ppm_topo(topoid)%t

         !----------------------------------------------------------------------
         !  scale the subs in this topology
         !----------------------------------------------------------------------
         DO j=1,topo%nsubs
             DO i=1,dim
#if __KIND == __DOUBLE_PRECISION
                topo%min_subd(i,j) = &
     &            (topo%min_subd(i,j) - origo(i))*scale(i) + origo(i)
                topo%max_subd(i,j) = &
     &           (topo%max_subd(i,j) - origo(i))*scale(i) + origo(i)
#else
                topo%min_subs(i,j) = &
     &            (topo%min_subs(i,j) - origo(i))*scale(i) + origo(i)
                topo%max_subs(i,j) = &
     &           (topo%max_subs(i,j) - origo(i))*scale(i) + origo(i)
#endif
		     ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !  scale the physical domain
         !----------------------------------------------------------------------
         DO i=1,dim
#if __KIND == __DOUBLE_PRECISION
             topo%min_physd(i) = &
     &         (topo%min_physd(i) - origo(i))*scale(i) + origo(i)
             topo%max_physd(i) = &
     &         (topo%max_physd(i) - origo(i))*scale(i) + origo(i)
#else
             topo%min_physs(i) = &
     &         (topo%min_physs(i) - origo(i))*scale(i) + origo(i)
             topo%max_physs(i) = &
     &         (topo%max_physs(i) - origo(i))*scale(i) + origo(i)
#endif
	    ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_scale_domain',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (topoid .GE. 0) THEN
            CALL ppm_check_topoid(topoid,valid,info)
            IF (.NOT. valid) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_scale_domain',  &
     &             'topoid is invalid!',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDIF
         !----------------------------------------------------------------------
         !  Check that the scale factor is greater than zero
         !----------------------------------------------------------------------
         DO k=1,ppm_dim
            IF (scale(k).LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_scale_domain',  &
     &         'the scale factor must > 0',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_scale_domain_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_scale_domain_d
#endif
