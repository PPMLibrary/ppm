      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_scale_domain
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine scales the computational domain ie. it 
      !                 scales the ppm_min_phys and ppm_max_phys and the subs 
      !                 associated with the topology. The scaling is given by:
      ! 
      !                    x_new = (x_old - x_origo)*x_scale + x_origo
      !
      !  Input        : topo_id      (I) user topology (not ppm internal!).
      !                                  topo_id => 0 scale the topoid
      !                                  topo_id = -1 scale all topologies
      !                 scale(:)     (F) scale factors 
      !                 origo(:)     (F) origin of the scaling
      !
      !  Input/output : 
      !
      !  Output       : info         (I) return status: zero on success
      !
      !  Remarks      : this routine is a user callable routine
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_scale_domain.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/09/04 18:34:54  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2004/08/31 14:15:38  hiebers
      !  removed typo in use module name
      !
      !  Revision 1.4  2004/08/31 13:29:58  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.3  2004/07/26 14:35:57  ivos
      !  Changed USE statements to new module structure.
      !
      !  Revision 1.2  2004/07/26 13:49:21  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.1  2004/07/19 14:21:56  walther
      !  First implementation (not tested).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_scale_domain_s(topo_id,scale,origo,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_scale_domain_d(topo_id,scale,origo,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_check_topoid
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
      REAL(MK), DIMENSION(:), INTENT(IN   ) :: scale,origo
      INTEGER               , INTENT(IN   ) :: topo_id ! user topoid
      INTEGER               , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
      REAL(MK)                          :: ox,oy,oz,sx,sy,sz
      INTEGER                           :: topoid ! ppm internal topoid
      INTEGER                           :: i,j,k,dim
      CHARACTER(LEN=ppm_char)           :: mesg
      LOGICAL                           :: valid
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
         IF (topo_id .GE. 0) THEN
            CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
            IF (.NOT. valid) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_scale_domain',  &
     &             'topo_id is invalid!',__LINE__,info)
               GOTO 9999
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
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the small vectors in registers for efficiency
      !-------------------------------------------------------------------------
      ox = origo(1)
      oy = origo(2)
      sx = scale(1)
      sy = scale(2)
      IF (ppm_dim.EQ.3) THEN
         oz = origo(3)
         sz = scale(3)
      ENDIF
 
      !-------------------------------------------------------------------------
      !  based on the topoid we have to scale all or one topology
      !-------------------------------------------------------------------------
      IF (topoid.LT.0) THEN
         !----------------------------------------------------------------------
         !  if the topology id is negative scale all topologies
         !----------------------------------------------------------------------
         DO k=1,ppm_max_topoid ! loop over all the topologies
            DO j=1,ppm_nsubs(k) ! loop over all subs in each k-topo
#if __KIND == __DOUBLE_PRECISION
               ppm_min_subd(1,j,k) = (ppm_min_subd(1,j,k) - ox)*sx + ox
               ppm_max_subd(1,j,k) = (ppm_max_subd(1,j,k) - ox)*sx + ox

               ppm_min_subd(2,j,k) = (ppm_min_subd(2,j,k) - oy)*sy + oy
               ppm_max_subd(2,j,k) = (ppm_max_subd(2,j,k) - oy)*sy + oy

               IF (ppm_dim.EQ.3) THEN
                  ppm_min_subd(3,j,k) = (ppm_min_subd(3,j,k) - oz)*sz + oz
                  ppm_max_subd(3,j,k) = (ppm_max_subd(3,j,k) - oz)*sz + oz
               ENDIF
#else
               ppm_min_subs(1,j,k) = (ppm_min_subs(1,j,k) - ox)*sx + ox
               ppm_max_subs(1,j,k) = (ppm_max_subs(1,j,k) - ox)*sx + ox

               ppm_min_subs(2,j,k) = (ppm_min_subs(2,j,k) - oy)*sy + oy
               ppm_max_subs(2,j,k) = (ppm_max_subs(2,j,k) - oy)*sy + oy

               IF (ppm_dim.EQ.3) THEN
                  ppm_min_subs(3,j,k) = (ppm_min_subs(3,j,k) - oz)*sz + oz
                  ppm_max_subs(3,j,k) = (ppm_max_subs(3,j,k) - oz)*sz + oz
               ENDIF
#endif
            ENDDO ! end loop of subs in k-th topo

            !-------------------------------------------------------------------
            !  scale the physical domain for/in each topology
            !-------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
            ppm_min_physd(1,k) = (ppm_min_physd(1,k) - ox)*sx + ox
            ppm_max_physd(1,k) = (ppm_max_physd(1,k) - ox)*sx + ox

            ppm_min_physd(2,k) = (ppm_min_physd(2,k) - oy)*sy + oy
            ppm_max_physd(2,k) = (ppm_max_physd(2,k) - oy)*sy + oy

            IF (ppm_dim.EQ.3) THEN
               ppm_min_physd(3,k) = (ppm_min_physd(3,k) - oz)*sz + oz
               ppm_max_physd(3,k) = (ppm_max_physd(3,k) - oz)*sz + oz
            ENDIF
#else
            ppm_min_physs(1,k) = (ppm_min_physs(1,k) - ox)*sx + ox
            ppm_max_physs(1,k) = (ppm_max_physs(1,k) - ox)*sx + ox

            ppm_min_physs(2,k) = (ppm_min_physs(2,k) - oy)*sy + oy
            ppm_max_physs(2,k) = (ppm_max_physs(2,k) - oy)*sy + oy

            IF (ppm_dim.EQ.3) THEN
               ppm_min_physs(3,k) = (ppm_min_physs(3,k) - oz)*sz + oz
               ppm_max_physs(3,k) = (ppm_max_physs(3,k) - oz)*sz + oz
            ENDIF
#endif
 
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  if the topology id (topo_id) is positive find the corresponding 
         !  internal ppm topoid (topoid)
         !----------------------------------------------------------------------
         topoid = ppm_internal_topoid(topo_id)

         !----------------------------------------------------------------------
         !  scale the subs in this topology
         !----------------------------------------------------------------------
         DO j=1,ppm_nsubs(topoid)
#if __KIND == __DOUBLE_PRECISION
            ppm_min_subd(1,j,topoid) = (ppm_min_subd(1,j,topoid) - ox)*sx + ox
            ppm_max_subd(1,j,topoid) = (ppm_max_subd(1,j,topoid) - ox)*sx + ox

            ppm_min_subd(2,j,topoid) = (ppm_min_subd(2,j,topoid) - oy)*sy + oy
            ppm_max_subd(2,j,topoid) = (ppm_max_subd(2,j,topoid) - oy)*sy + oy
            IF (ppm_dim.EQ.3) THEN
               ppm_min_subd(3,j,topoid) = (ppm_min_subd(3,j,topoid) - oz)*sz +oz
               ppm_max_subd(3,j,topoid) = (ppm_max_subd(3,j,topoid) - oz)*sz +oz
            ENDIF 
#else
            ppm_min_subs(1,j,topoid) = (ppm_min_subs(1,j,topoid) - ox)*sx + ox
            ppm_max_subs(1,j,topoid) = (ppm_max_subs(1,j,topoid) - ox)*sx + ox

            ppm_min_subs(2,j,topoid) = (ppm_min_subs(2,j,topoid) - oy)*sy + oy
            ppm_max_subs(2,j,topoid) = (ppm_max_subs(2,j,topoid) - oy)*sy + oy
            IF (ppm_dim.EQ.3) THEN
               ppm_min_subs(3,j,topoid) = (ppm_min_subs(3,j,topoid) - oz)*sz +oz
               ppm_max_subs(3,j,topoid) = (ppm_max_subs(3,j,topoid) - oz)*sz +oz
            ENDIF
#endif
         ENDDO

         !----------------------------------------------------------------------
         !  scale the physical domain
         !----------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
         ppm_min_physd(1,topoid) = (ppm_min_physd(1,topoid) - ox)*sx + ox
         ppm_max_physd(1,topoid) = (ppm_max_physd(1,topoid) - ox)*sx + ox

         ppm_min_physd(2,topoid) = (ppm_min_physd(2,topoid) - oy)*sy + oy
         ppm_max_physd(2,topoid) = (ppm_max_physd(2,topoid) - oy)*sy + oy
         IF (ppm_dim.EQ.3) THEN
            ppm_min_physd(3,topoid) = (ppm_min_physd(3,topoid) - oz)*sz + oz
            ppm_max_physd(3,topoid) = (ppm_max_physd(3,topoid) - oz)*sz + oz
         ENDIF
#else
         ppm_min_physs(1,topoid) = (ppm_min_physs(1,topoid) - ox)*sx + ox
         ppm_max_physs(1,topoid) = (ppm_max_physs(1,topoid) - ox)*sx + ox

         ppm_min_physs(2,topoid) = (ppm_min_physs(2,topoid) - oy)*sy + oy
         ppm_max_physs(2,topoid) = (ppm_max_physs(2,topoid) - oy)*sy + oy
         IF (ppm_dim.EQ.3) THEN
            ppm_min_physs(3,topoid) = (ppm_min_physs(3,topoid) - oz)*sz + oz
            ppm_max_physs(3,topoid) = (ppm_max_physs(3,topoid) - oz)*sz + oz
         ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_scale_domain',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_scale_domain_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_scale_domain_d
#endif
