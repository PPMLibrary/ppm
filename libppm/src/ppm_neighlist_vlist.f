      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_neighlist_vlist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Create Verlet lists for all particles of this processor.
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 Np         (I) number of particles on this proc.
      !                 cutoff     (F) cutoff radius for PP interactions
      !                 skin       (F) Verlet list skin layer thikness.
      !                 lsymm      (L) T to use symmetry, F to not do so.
      !                 pidx(:)    (I) indices of those particles that are
      !                                to be included in the list
      !                                (OPTIONAL). By default all particles
      !                                are taken. If given, particles
      !                                indices in Verlet lists are relative
      !                                to xp(:,pidx) and not xp(:,:)!!
      !                 lstore     (L) OTIONAL Set this to .TRUE. to store 
      !                                (and return) the Verlet lists in
      !                                vlist. If this is false, only
      !                                nvlist is determined and
      !                                returned. Default is .TRUE.
      !
      !  Input/output : 
      !
      !  Output       : info       (I) return status. =0 if no error.
      !                 vlist(:,:) (I) Verlet list. First index: particles
      !                                with which particle ip interacts.
      !                                Second index: ip.
      !                                The second index only runs up to the
      !                                largest ip with non-zero nvlist.
      !                                This is to save memory since the
      !                                last particles are the ghosts and
      !                                they do not have a verlet list.
      !                                This is only allocated and
      !                                returned if lstore is .TRUE.
      !                 nvlist(:)  (I) Number of particles with which ip
      !                                has to interact. Index: ip.
      !
      !  Remarks      : The list needs to be rebuilt as soon as a particle
      !                 has moved a distance larger than skin. It is the
      !                 USERs responsibility to detect when this is the
      !                 case and CALL this routine again.
      !
      !                 vlist and nvlist are allocated in this routine. The
      !                 user just needs to pass pointers.
      !
      !                 the two cases for lsymm have their own duplicated
      !                 loops since the lsymm=F case does not vectorize.
      !                 lsymm=T (using symmetry) however does.
      !
      !                 The __VECTOR case was tested and found to vectorize
      !                 on the NEC SX-5 even without compiler directives.
      !
      !                 Requires (almost) two repetitions of the main loops.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_neighlist_vlist.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.26  2006/07/04 15:58:59  ivos
      !  Added comments and changed alloc_fit to alloc_grow.
      !
      !  Revision 1.25  2006/07/04 15:24:51  ivos
      !  Added vectorization remark to header comment.
      !
      !  Revision 1.24  2004/11/11 15:26:20  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.23  2004/11/09 11:58:09  ivos
      !  Added OPTIONAL argument lstore.
      !
      !  Revision 1.22  2004/11/04 16:28:54  ivos
      !  Optimized loops for the case of empty cells.
      !
      !  Revision 1.21  2004/10/01 16:09:11  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.20  2004/07/26 15:38:50  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.19  2004/07/26 11:59:39  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.18  2004/07/26 07:42:49  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.17  2004/07/16 14:47:21  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.16  2004/07/15 09:02:41  ivos
      !  bugfix: initialized kk to 0 instead of nvlist(ip).
      !
      !  Revision 1.15  2004/07/14 09:31:04  ivos
      !  Added VECTOR stuff and made all inner loops vectorize.
      !
      !  Revision 1.14  2004/07/14 07:11:25  ivos
      !  bugfix: double-interactions in the box itself added for the non-
      !  symmetric case.
      !
      !  Revision 1.13  2004/06/15 08:43:30  ivos
      !  bugfix: Nm is now a 2d array since it depends on the subdomain id.
      !  n1 and n2 are now computed inside the loop. Also added SXF90 directives
      !  for vectorization.
      !
      !  Revision 1.12  2004/06/03 16:08:34  ivos
      !  Fixed over-tight argument check skin.LE.0 to skin.LT.0 and added
      !  dealloc of inp and jnp at the end.
      !
      !  Revision 1.11  2004/02/24 11:39:05  ivos
      !  bugfix: changed CALL to clist to new argument syntax.
      !
      !  Revision 1.10  2004/02/24 11:36:39  ivos
      !  Changed imode (INTEGER symmetry flag) argument to lsymm (LOGICAL) in
      !  order to have the same interface as the ghost routines.
      !
      !  Revision 1.9  2004/02/04 17:20:48  ivos
      !  renamed TYPE ptr_to_clist to ppm_type_ptr_to_clist.
      !
      !  Revision 1.8  2004/01/27 14:50:36  ivos
      !  Attempt to make routine vectorize: (1) removed calls to ppm_alloc from
      !  inner loops. Instead, upper bound to list length is estimated first and
      !  then the lists are pruned just before exit. The total number of 
      !  particles in all interacting cells is taken as an upper bound. 
      !  (2) removed inner-most
      !  loop over ppm_dim and replaced with explicit dx,dy,dz. (3) Added
      !  OPTIONAL argument to specify a subset list of the particles.
      !  PRELIMINARY COMMIT. NONE OF THIS HAS BEEN TESTED YET.
      !
      !  Revision 1.7  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.6  2004/01/22 14:55:45  ivos
      !  Added checks after allocs.
      !
      !  Revision 1.5  2004/01/22 13:27:56  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.4  2004/01/12 12:38:37  ivos
      !  Uses ppm_neighlist_clist_destroy now to properly dealloc cell list.
      !
      !  Revision 1.3  2004/01/09 16:29:19  ivos
      !  Changed to new interaction handling and added flag to switch symmetry
      !  on/off. duplicated the loops and moved the IF(lsymm) outside to
      !  allow the symmetric case to vectorize.
      !
      !  Revision 1.2  2004/01/08 17:53:34  ivos
      !  Cosmetics.
      !
      !  Revision 1.1  2004/01/08 17:52:14  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_neighlist_vlist_s(xp,Np,cutoff,skin,lsymm,vlist, &
     &               nvlist,info,pidx,lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_neighlist_vlist_d(xp,Np,cutoff,skin,lsymm,vlist, &
     &               nvlist,info,pidx,lstore)
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_neighlist
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_neighlist_MkNeighIdx
      USE ppm_module_neighlist_clist
      USE ppm_module_clist_destroy
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      INTEGER                 , INTENT(IN   ) :: Np
      REAL(MK)                , INTENT(IN   ) :: cutoff,skin
      LOGICAL                 , INTENT(IN   ) :: lsymm
      INTEGER                 , INTENT(  OUT) :: info
      INTEGER, DIMENSION(:,:) , POINTER       :: vlist
      INTEGER, DIMENSION(  :) , POINTER       :: nvlist
      INTEGER, DIMENSION(  :) , OPTIONAL      :: pidx
      LOGICAL, INTENT(IN)     , OPTIONAL      :: lstore
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                   :: t0
      ! effective number of particles
      INTEGER                                    :: Npdx
      ! counters
      INTEGER                                    :: i,idom,ibox,jbox,nbox
      INTEGER                                    :: ipart,jpart,ip,jp,maxvlen
      INTEGER                                    :: cbox,iinter,j,k
      INTEGER                                    :: ii,jj,kk
      ! coordinate difference
      REAL(MK)                                   :: dx,dy,dz
      ! inter particle distance
      REAL(MK)                                   :: dij
      ! cutoff squared
      REAL(MK)                                   :: cut2
      ! start and end particle in a box
      INTEGER                                    :: istart,iend,jstart,jend
      ! box size for helper cell list
      REAL(MK), DIMENSION(3)                     :: bsize
      ! cell neighbor lists
      INTEGER, DIMENSION(:,:), POINTER           :: inp,jnp
      ! number of interactions for each cell
      INTEGER                                    :: nnp
      ! for allocate
      INTEGER, DIMENSION(2)                      :: lda
      INTEGER                                    :: iopt
      ! number of cells in all directions
      INTEGER, DIMENSION(:,:), POINTER           :: Nm
      ! cell offsets for box index
      INTEGER                                    :: n1,n2,nz,lb
      CHARACTER(LEN=ppm_char)                    :: mesg
      ! store vlist?
      LOGICAL                                    :: lst
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_neighlist_vlist',t0,info)
      !-------------------------------------------------------------------------
      !  If the user gave an explicit list of particles to be included, use
      !  the size of this list as the effective number of particles. Use
      !  Np otherwise.
      !-------------------------------------------------------------------------
      IF (PRESENT(pidx)) THEN
          IF (Np .GT. SIZE(pidx,1)) Npdx = SIZE(pidx,1)
      ELSE
          Npdx = Np
      ENDIF
      !-------------------------------------------------------------------------
      !  Do we need to store the Verlet lists or just determine their lengths?
      !-------------------------------------------------------------------------
      IF (PRESENT(lstore)) THEN
          lst = lstore
      ELSE
          lst = .TRUE.
      ENDIF
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_neighlist_vlist',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (cutoff .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'cutoff must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (skin .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'skin must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Npdx .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'Np must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Boxes need to be cutoff+skin in all directions !
      !-------------------------------------------------------------------------
      DO i=1,ppm_dim
          bsize(i) = cutoff + skin
      ENDDO
      cut2 = bsize(1)*bsize(1)
      !-------------------------------------------------------------------------
      !  Generate cell lists 
      !-------------------------------------------------------------------------
      IF (PRESENT(pidx)) THEN
          CALL ppm_neighlist_clist(xp(:,pidx),Npdx,bsize,lsymm,clist,Nm,info)
      ELSE
          CALL ppm_neighlist_clist(xp,Npdx,bsize,lsymm,clist,Nm,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_vlist',   &
     &          'Building cell lists failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Generate cell neighbor lists 
      !-------------------------------------------------------------------------
      CALL ppm_neighlist_MkNeighIdx(lsymm,inp,jnp,nnp,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_vlist',  &
     &        'neighlist_MkNeighIdx failed',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate memory for Verlet list lengths
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
      lda(1) = Npdx
      CALL ppm_alloc(nvlist,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_neighlist_vlist',  &
     &         'Verlet list sizes NVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      nvlist(1:Npdx) = 0

      !-------------------------------------------------------------------------
      !  Lower box bound depends on symmetry use
      !-------------------------------------------------------------------------
      IF (lsymm) THEN 
          lb = 0
      ELSE
          lb = 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine size of Verlet lists
      !-------------------------------------------------------------------------
      DO idom=1,ppm_nsublist(ppm_topoid)
          n1  = Nm(1,idom)
          n2  = Nm(1,idom)*Nm(2,idom)
          nz  = Nm(3,idom)
          IF (ppm_dim .EQ. 2) THEN
              n2 = 0
              nz = 2
          ENDIF 
          ! loop over all REAL cells (the -2 at the end does this)
          DO k=lb,nz-2
              DO j=lb,Nm(2,idom)-2
                  DO i=lb,Nm(1,idom)-2
                      ! index of the center box
                      cbox = i + 1 + n1*j + n2*k
                      ! loop over all box-box interactions
                      DO iinter=1,nnp
                          ! determine box indices for this interaction
                          ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                           n2*inp(3,iinter))
                          jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                           n2*jnp(3,iinter))
                          !-----------------------------------------------------
                          !  Read indices and check if empty
                          !-----------------------------------------------------
                          istart = clist(idom)%lhbx(ibox)
                          iend   = clist(idom)%lhbx(ibox+1)-1
                          IF (iend .LT. istart) CYCLE
                          !-----------------------------------------------------
                          !  Within the box itself use symmetry and avoid 
                          !  adding the particle itself to its own list
                          !-----------------------------------------------------
                          IF (ibox .EQ. jbox) THEN
                              DO ipart=istart,iend
                                  ip = clist(idom)%lpdx(ipart)
                                  IF (lsymm) THEN
                                      DO jpart=(ipart+1),iend
                                          jp = clist(idom)%lpdx(jpart)
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to list of
                                          ! particle ip
                                          nvlist(ip) = nvlist(ip) + 1
                                      ENDDO
                                  ELSE
#ifdef __VECTOR
                                      DO jpart=istart,iend
                                          jp = clist(idom)%lpdx(jpart)
                                          IF (jp .EQ. ip) CYCLE
#else
                                      DO jpart=(ipart+1),iend
                                          jp = clist(idom)%lpdx(jpart)
#endif
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to list of
                                          ! particle ip
                                          nvlist(ip) = nvlist(ip) + 1
#ifndef __VECTOR
                                          nvlist(jp) = nvlist(jp) + 1
#endif
                                      ENDDO
                                  ENDIF
                              ENDDO
                          !-----------------------------------------------------
                          !  For the other boxes check all particles
                          !-----------------------------------------------------
                          ELSE
                              ! get pointers to first and last particle 
                              jstart = clist(idom)%lhbx(jbox)
                              jend   = clist(idom)%lhbx(jbox+1)-1
                              ! skip this iinter if empty
                              If (jend .LT. jstart) CYCLE
                              ! loop over all particles inside this cell
                              DO ipart=istart,iend
                                  ip = clist(idom)%lpdx(ipart)
                                  ! check against all particles 
                                  ! in the other cell
                                  DO jpart=jstart,jend
                                      jp = clist(idom)%lpdx(jpart)
                                      ! translate to real particle
                                      ! index if needed
                                      IF (PRESENT(pidx)) THEN
                                          ii = pidx(ip)
                                          jj = pidx(jp)
                                      ELSE
                                          ii = ip
                                          jj = jp
                                      ENDIF
                                      dx = xp(1,ii) - xp(1,jj)
                                      dy = xp(2,ii) - xp(2,jj)
                                      IF (ppm_dim .GT. 2) THEN
                                          dz = xp(3,ii) - xp(3,jj)
                                          dij= (dx*dx)+(dy*dy)+(dz*dz)
                                      ELSE
                                          dz = 0.0_MK
                                          dij= (dx*dx)+(dy*dy)
                                      ENDIF
                                      IF (dij .GT. cut2) CYCLE
                                      ! add particle jp to Verlet 
                                      ! list of particle ip
                                      nvlist(ip) = nvlist(ip) + 1
                                  ENDDO
                              ENDDO
                          ENDIF       ! ibox .EQ. jbox
                      ENDDO           ! iinter
                  ENDDO               ! i
              ENDDO                   ! j
          ENDDO                       ! k
      ENDDO                           ! idom

      !-------------------------------------------------------------------------
      !  Maximum Verlet list length
      !-------------------------------------------------------------------------
      maxvlen = MAXVAL(nvlist)
      IF (ppm_debug .GT. 0) THEN
          WRITE(mesg,'(A,I8)') 'Maximum length of Verlet lists: ',maxvlen
          CALL ppm_write(ppm_rank,'ppm_neighlist_vlist',mesg,info)
      ENDIF
      !-------------------------------------------------------------------------
      !  Only build and store the lists if needed
      !-------------------------------------------------------------------------
      IF (lst) THEN
          !---------------------------------------------------------------------
          !  Highest particle index with non-zero nvlist
          !---------------------------------------------------------------------
          ii = Npdx
          DO ipart=Npdx,1,-1
              IF (nvlist(ipart) .GT. 0) THEN
                  ii = ipart
                  EXIT
              ENDIF
          ENDDO
          !---------------------------------------------------------------------
          !  Allocate memory for Verlet lists
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_grow
          lda(1) = maxvlen
          lda(2) = ii
          CALL ppm_alloc(vlist,lda,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_vlist',  &
         &         'Verlet list VLIST',__LINE__,info)
              GOTO 9999
          ENDIF
          nvlist(1:Npdx) = 0
          !---------------------------------------------------------------------
          !  BUILD VERLET LISTS 
          !---------------------------------------------------------------------
          DO idom=1,ppm_nsublist(ppm_topoid)
              n1  = Nm(1,idom)
              n2  = Nm(1,idom)*Nm(2,idom)
              nz  = Nm(3,idom)
              IF (ppm_dim .EQ. 2) THEN
                  n2 = 0
                  nz = 2
              ENDIF 
              ! get number of cells in this subdomain
              nbox = SIZE(clist(idom)%lhbx,1)-1
              ! loop over all REAL cells (the -2 at the end does this)
              DO k=lb,nz-2
                  DO j=lb,Nm(2,idom)-2
                      DO i=lb,Nm(1,idom)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                           n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                           n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if empty
                              !-------------------------------------------------
                              istart = clist(idom)%lhbx(ibox)
                              iend   = clist(idom)%lhbx(ibox+1)-1
                              IF (iend .LT. istart) CYCLE
                              !-------------------------------------------------
                              !  Within the box itself use symmetry and avoid 
                              !  adding the particle itself to its own list
                              !-------------------------------------------------
                              IF (ibox .EQ. jbox) THEN
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      kk = nvlist(ip)
                                      IF (lsymm) THEN
                                          DO jpart=(ipart+1),iend
                                              jp = clist(idom)%lpdx(jpart)
                                              ! translate to real particle
                                              ! index if needed
                                              IF (PRESENT(pidx)) THEN
                                                  ii = pidx(ip)
                                                  jj = pidx(jp)
                                              ELSE
                                                  ii = ip
                                                  jj = jp
                                              ENDIF
                                              dx = xp(1,ii) - xp(1,jj)
                                              dy = xp(2,ii) - xp(2,jj)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz = xp(3,ii) - xp(3,jj)
                                                  dij= (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij= (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cut2) CYCLE
                                              ! add particle jp to list of
                                              ! particle ip
                                              kk = kk + 1
                                              vlist(kk,ip) = jp
                                          ENDDO
                                      ELSE
#ifdef __VECTOR
                                          DO jpart=istart,iend
                                              jp = clist(idom)%lpdx(jpart)
                                              IF (jp .EQ. ip) CYCLE
#else
                                          DO jpart=(ipart+1),iend
                                              jp = clist(idom)%lpdx(jpart)
#endif
                                              ! translate to real particle
                                              ! index if needed
                                              IF (PRESENT(pidx)) THEN
                                                  ii = pidx(ip)
                                                  jj = pidx(jp)
                                              ELSE
                                                  ii = ip
                                                  jj = jp
                                              ENDIF
                                              dx = xp(1,ii) - xp(1,jj)
                                              dy = xp(2,ii) - xp(2,jj)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz = xp(3,ii) - xp(3,jj)
                                                  dij= (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij= (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cut2) CYCLE
                                              kk = kk + 1
                                              ! add particle jp to list of
                                              ! particle ip
                                              vlist(kk,ip) = jp
#ifndef __VECTOR
                                              ! add particle ip to list of
                                              ! particle jp
                                              nvlist(jp) = nvlist(jp) + 1
                                              vlist(nvlist(jp),jp) = ip
#endif
                                          ENDDO
                                      ENDIF
                                      nvlist(ip) = kk
                                  ENDDO
                              !-------------------------------------------------
                              !  For the other boxes check all particles
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = clist(idom)%lhbx(jbox)
                                  jend   = clist(idom)%lhbx(jbox+1)-1
                                  ! skip this iinter if empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      kk = nvlist(ip)
                                      ! check against all particles 
                                      ! in the other cell
                                      DO jpart=jstart,jend
                                          jp = clist(idom)%lpdx(jpart)
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to Verlet 
                                          ! list of particle ip
                                          kk = kk + 1
                                          vlist(kk,ip) = jp
                                      ENDDO
                                      nvlist(ip) = kk
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom
      ENDIF       ! lstore

      !-------------------------------------------------------------------------
      !  Free work cell lists
      !-------------------------------------------------------------------------
      CALL ppm_clist_destroy(clist,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Failed to destroy cell lists',__LINE__,info)
      ENDIF
      iopt = ppm_param_dealloc
      CALL ppm_alloc(inp,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Box interaction index INP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(jnp,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Box interaction index JNP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(Nm,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Numbers of cells NM',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_neighlist_vlist',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_neighlist_vlist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_neighlist_vlist_d
#endif
