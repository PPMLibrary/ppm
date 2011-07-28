      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_alloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine (re-)allocates a topology object 
      !                 and all its members and marks it as undefined.
      !
      !  Input        : nsubs      (I) Total number of subs on all proc
      !                 nsublist   (I) Local number of subs on this proc
      !                 maxneigh   (I) maximum number of neighbors of any
      !                                sub on this processor
      !                 prec       (I) precision for storage. One of:
      !                                    ppm_kind_single
      !                                    ppm_kind_double
      !
      !  Input/output : topo  (ppm_t_topo) Topology structure to be
      !                                allocated
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      : 
      !                
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_alloc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:08  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_topo_alloc(topo,nsubs,nsublist,maxneigh,prec,info)

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
      TYPE(ppm_t_topo)        , INTENT(INOUT) :: topo
      INTEGER                 , INTENT(IN   ) :: nsubs,nsublist
      INTEGER                 , INTENT(IN   ) :: maxneigh,prec
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3) :: ldc
      INTEGER                :: iopt,nsubmax
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_alloc',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_alloc',  &
     &            'nsubs must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsublist .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_alloc',  &
     &            'nsublist must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsubs .LE. nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_alloc',  &
     &            'total number of subs is smaller than local',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (maxneigh .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_alloc',  &
     &            'maxneigh must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  The MAX of naublist and 1 is needed to avoid allocation failures
      !  if a processor has 0 subs
      !-------------------------------------------------------------------------
      nsubmax = MAX(nsublist,1)

      !-------------------------------------------------------------------------
      !  Allocate memory for the domain
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      IF (prec .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(topo%min_physs,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'min extent of domain TOPO%MIN_PHYSS',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(topo%max_physs,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'max extent of domain TOPO%MAX_PHYSS',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ASSOCIATED(topo%min_physd)) DEALLOCATE(topo%min_physd,STAT=info)
          NULLIFY(topo%min_physd)
          IF (ASSOCIATED(topo%max_physd)) DEALLOCATE(topo%max_physd,STAT=info)
          NULLIFY(topo%max_physd)
      ELSE
          CALL ppm_alloc(topo%min_physd,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'min extent of domain TOPO%MIN_PHYSD',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(topo%max_physd,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'max extent of domain TOPO%MAX_PHYSD',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ASSOCIATED(topo%min_physs)) DEALLOCATE(topo%min_physs,STAT=info)
          NULLIFY(topo%min_physs)
          IF (ASSOCIATED(topo%max_physs)) DEALLOCATE(topo%max_physs,STAT=info)
          NULLIFY(topo%max_physs)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the subdomains 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      IF (prec .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(topo%min_subs,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'min extent of subs TOPO%MIN_SUBS',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(topo%max_subs,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'max extent of subs TOPO%MAX_SUBS',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ASSOCIATED(topo%min_subd)) DEALLOCATE(topo%min_subd,STAT=info)
          NULLIFY(topo%min_subd)
          IF (ASSOCIATED(topo%max_subd)) DEALLOCATE(topo%max_subd,STAT=info)
          NULLIFY(topo%max_subd)
      ELSE
          CALL ppm_alloc(topo%min_subd,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'min extent of subs TOPO%MIN_SUBD',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(topo%max_subd,ldc,iopt,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &            'max extent of subs TOPO%MAX_SUBD',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ASSOCIATED(topo%min_subs)) DEALLOCATE(topo%min_subs,STAT=info)
          NULLIFY(topo%min_subs)
          IF (ASSOCIATED(topo%max_subs)) DEALLOCATE(topo%max_subs,STAT=info)
          NULLIFY(topo%max_subs)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the sub-to-proc mapping
      !-------------------------------------------------------------------------
      ldc(1) = nsubs
      CALL ppm_alloc(topo%subs2proc,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'global subs to proc map TOPO%SUBS2PROC',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for list of local subs 
      !-------------------------------------------------------------------------
      ldc(1) = nsubmax
      CALL ppm_alloc(topo%isublist,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'local sub list TOPO%ISUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the boundary conditions
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 2*ppm_dim
      CALL ppm_alloc(topo%bcdef,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'boundary conditions TOPO%BCDEF',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Allocate memory for the subdomain costs
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs
      IF (prec .EQ. ppm_kind_single) THEN
          IF (ASSOCIATED(topo%sub_costd)) DEALLOCATE(topo%sub_costd,STAT=info)
          NULLIFY(topo%sub_costd)
          CALL ppm_alloc(topo%sub_costs,ldc,iopt,info)
      ELSE
          IF (ASSOCIATED(topo%sub_costs)) DEALLOCATE(topo%sub_costs,STAT=info)
          NULLIFY(topo%sub_costs)
          CALL ppm_alloc(topo%sub_costd,ldc,iopt,info)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'subdomain costs TOPO%SUB_COST',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Allocate memory for the numbers of neighbors of each local sub on
      !  the current processor i.e, a subset of the nneigh(1:nsubs) list.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubmax
      CALL ppm_alloc(topo%nneighsubs,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'number of neighbors of local subs TOPO%NNEIGHSUBS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global IDs of all the neighbors.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = maxneigh
      ldc(2) = nsubmax
      CALL ppm_alloc(topo%ineighsubs,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'neighbors of local subs TOPO%INEIGHSUBS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for BC for the subs on the current processor
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 2*ppm_dim
      ldc(2) = nsubmax
      CALL ppm_alloc(topo%subs_bc,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'BCs for subs on local processor TOPO%SUBS_BC',__LINE__,info)
          GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate memory for the list of neighboring processors.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 26    ! just guessing :-)
      CALL ppm_alloc(topo%ineighproc,ldc,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_alloc',     &
     &        'list of neighbor processors TOPO%INEIGHPROC',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Nullify the communication sequence. It will be allocated when
      !  first used.
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(topo%icommseq)) DEALLOCATE(topo%icommseq,STAT=info)
      NULLIFY(topo%icommseq)
      topo%isoptimized = .FALSE.

      !-------------------------------------------------------------------------
      !  By default there are no meshes defined on a topology
      !-------------------------------------------------------------------------
      topo%max_meshid = 0
      iopt = ppm_param_dealloc
      CALL ppm_mesh_alloc_equi(topo%mesh,ldc,iopt,info)
      NULLIFY(topo%mesh)

      !-------------------------------------------------------------------------
      !  Mark this topology as not defined
      !-------------------------------------------------------------------------
      topo%isdefined = .FALSE.

      !-------------------------------------------------------------------------
      !  Set the precision for this topology
      !-------------------------------------------------------------------------
      topo%prec = prec

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_alloc',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_alloc
