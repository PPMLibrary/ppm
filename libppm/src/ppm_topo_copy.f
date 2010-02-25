      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_copy
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine copies the member structures from
      !                 one topology object to another one. The previous
      !                 contents of outtopo will be destroyed. This also
      !                 allocates outtopo.
      !
      !  Input        : intopo  (ppm_t_topo) Input topology
      !
      !  Input/output : 
      !
      !  Output       : outtopo (ppm_t_topo) Output topology. Will contain
      !                                the same data as intopo upon return.
      !                 info        (I) return status. 0 upon success.
      !
      !  Remarks      : 
      !                
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_copy.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:08  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_topo_copy(intopo,outtopo,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo)        , INTENT(IN   ) :: intopo
      TYPE(ppm_t_topo)        , INTENT(  OUT) :: outtopo
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)  :: t0
      INTEGER                :: nsubs,nsublist,maxneigh,prec,i,j,iopt
      INTEGER, DIMENSION(1)  :: ldc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_copy',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
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
      CALL ppm_topo_alloc(outtopo,nsubs,nsublist,maxneigh,prec,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_copy',     &
     &        'result topology OUTTOPO',__LINE__,info)
          GOTO 9999
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
      outtopo%max_meshid  = intopo%max_meshid
      
      !-------------------------------------------------------------------------
      !  Copy the boundary conditions
      !-------------------------------------------------------------------------
      DO i=1,2*ppm_dim
          outtopo%bcdef(i) = intopo%bcdef(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy the domain
      !-------------------------------------------------------------------------
      IF (prec .EQ. ppm_kind_single) THEN
          IF (ppm_dim .EQ. 2) THEN
              outtopo%min_physs(1) = intopo%min_physs(1)
              outtopo%min_physs(2) = intopo%min_physs(2)
              outtopo%max_physs(1) = intopo%max_physs(1)
              outtopo%max_physs(2) = intopo%max_physs(2)
          ELSE
              outtopo%min_physs(1) = intopo%min_physs(1)
              outtopo%min_physs(2) = intopo%min_physs(2)
              outtopo%min_physs(3) = intopo%min_physs(3)
              outtopo%max_physs(1) = intopo%max_physs(1)
              outtopo%max_physs(2) = intopo%max_physs(2)
              outtopo%max_physs(3) = intopo%max_physs(3)
          ENDIF
      ELSE
          IF (ppm_dim .EQ. 2) THEN
              outtopo%min_physd(1) = intopo%min_physd(1)
              outtopo%min_physd(2) = intopo%min_physd(2)
              outtopo%max_physd(1) = intopo%max_physd(1)
              outtopo%max_physd(2) = intopo%max_physd(2)
          ELSE
              outtopo%min_physd(1) = intopo%min_physd(1)
              outtopo%min_physd(2) = intopo%min_physd(2)
              outtopo%min_physd(3) = intopo%min_physd(3)
              outtopo%max_physd(1) = intopo%max_physd(1)
              outtopo%max_physd(2) = intopo%max_physd(2)
              outtopo%max_physd(3) = intopo%max_physd(3)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the subdomains
      !-------------------------------------------------------------------------
      IF (prec .EQ. ppm_kind_single) THEN
          IF (ppm_dim .EQ. 2) THEN
              DO i=1,nsubs
                  outtopo%min_subs(1,i) = intopo%min_subs(1,i)
                  outtopo%min_subs(2,i) = intopo%min_subs(2,i)
                  outtopo%max_subs(1,i) = intopo%max_subs(1,i)
                  outtopo%max_subs(2,i) = intopo%max_subs(2,i)
              ENDDO
          ELSE
              DO i=1,nsubs
                  outtopo%min_subs(1,i) = intopo%min_subs(1,i)
                  outtopo%min_subs(2,i) = intopo%min_subs(2,i)
                  outtopo%min_subs(3,i) = intopo%min_subs(3,i)
                  outtopo%max_subs(1,i) = intopo%max_subs(1,i)
                  outtopo%max_subs(2,i) = intopo%max_subs(2,i)
                  outtopo%max_subs(3,i) = intopo%max_subs(3,i)
              ENDDO
          ENDIF
      ELSE
          IF (ppm_dim .EQ. 2) THEN
              DO i=1,nsubs
                  outtopo%min_subd(1,i) = intopo%min_subd(1,i)
                  outtopo%min_subd(2,i) = intopo%min_subd(2,i)
                  outtopo%max_subd(1,i) = intopo%max_subd(1,i)
                  outtopo%max_subd(2,i) = intopo%max_subd(2,i)
              ENDDO
          ELSE
              DO i=1,nsubs
                  outtopo%min_subd(1,i) = intopo%min_subd(1,i)
                  outtopo%min_subd(2,i) = intopo%min_subd(2,i)
                  outtopo%min_subd(3,i) = intopo%min_subd(3,i)
                  outtopo%max_subd(1,i) = intopo%max_subd(1,i)
                  outtopo%max_subd(2,i) = intopo%max_subd(2,i)
                  outtopo%max_subd(3,i) = intopo%max_subd(3,i)
              ENDDO
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the sub costs
      !-------------------------------------------------------------------------
      IF (prec .EQ. ppm_kind_single) THEN
          DO i=1,nsubs
              outtopo%sub_costs(i) = intopo%sub_costs(i)
          ENDDO
      ELSE
          DO i=1,nsubs
              outtopo%sub_costd(i) = intopo%sub_costd(i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the sub to proc assignment
      !-------------------------------------------------------------------------
      DO i=1,nsubs
          outtopo%subs2proc(i) = intopo%subs2proc(i)
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
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_copy',     &
     &        'list of neighboring proc OUTTOPO%INEIGHPROC',__LINE__,info)
          GOTO 9999
      ENDIF
      DO i=1,intopo%nneighproc
          outtopo%ineighproc(i) = intopo%ineighproc(i)
      ENDDO

      IF (intopo%isoptimized) THEN
          iopt   = ppm_param_alloc_fit
          ldc(1) = intopo%ncommseq
          CALL ppm_alloc(outtopo%icommseq,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_copy',     &
         &        'communication sequence OUTTOPO%ICOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,intopo%ncommseq
              outtopo%icommseq(i) = intopo%icommseq(i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the mesh definitions
      !-------------------------------------------------------------------------
      IF (intopo%max_meshid .GT. 0) THEN
          ALLOCATE(outtopo%mesh(intopo%max_meshid),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_copy',     &
         &        'mesh definitions OUTTOPO%MESH',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,intopo%max_meshid
              outtopo%mesh(i) = intopo%mesh(i)
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
      CALL substop('ppm_topo_copy',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_copy
