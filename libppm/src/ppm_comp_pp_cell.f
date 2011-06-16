      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_comp_pp_cell
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine computes kernel interactions by
      !                 direct particle-particle interactions using cell
      !                 lists.
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 Np         (I) number of particles on this proc.
      !                 pdata(:,:) (O) particle strengths used for interaction. 
      !                                Overloaded types:
      !                                single,double,single complex and
      !                                double complex.
      !                 lda        (I) lading dimension of pdata.
      !                 lsymm      (L) Do symmetric PP interactions?
      !                                   .TRUE.  : Yes
      !                                   .FALSE. : No
      !                 kernel     (O) kernel to be used for PP
      !                                interactions. To use ppm-internal
      !                                kernels, specify one of:
      !                                   ppm_param_kerel_laplace2d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 2D)
      !                                   ppm_param_kerel_laplace3d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 3D)
      !                                To use your own kernel function,
      !                                pass the function pointer here. Your
      !                                function should take one argument
      !                                and return one value. The third
      !                                possibility is to pass a lookup
      !                                table with tabulated kernel values.
      !                                Such a table can be created using
      !                                ppm_comp_pp_mk_table.
      !                 kpar(:)    (O) Kernel parameters. See documentation
      !                                or ppm_comp_pp_kernels.inc for
      !                                description. Type can be single,
      !                                double, single complex or double
      !                                complex. When using a lookup table,
      !                                pass dxtableinv (the inverse of the
      !                                table spacing) as a scalar here.
      !                 nsublist   (I) number of subdomains on the local
      !                                processor
      !                 clist(:)       Cell list as a list of ptr_to_clist.
      !                                particle index list of isub: 
      !                                        clist(isub)%lpdx(:)
      !                                pointer to first particle in ibox of
      !                                isub: 
      !                                        clist(isub)%lhbx(ibox)
      !                 Nm(:,:)    (I) number of cells in x,y,(z)
      !                                direction (including the ghosts cells) 
      !                                in each subdomain. 1st index:
      !                                direction. second index: subid.
      !                 cutoff2    (F) Squared PP interaction cutoff.
      !                                Should be .LE. cell size.
      !
      !  Input/output : 
      !
      !  Output       : dpd(:,:)   (O) Change of particle data (pdata) due to 
      !                                interaction. This is not initialized
      !                                by this routine! Overloaded types:
      !                                single,double,single complex,double
      !                                complex.
      !                 info       (I) return status. =0 if no error.
      !
      !  Remarks      : Loops over lda have been explicitly unrolled for
      !                 the cases of lda.EQ.1 and lda.EQ.2 to allow
      !                 vectorization in these cases.
      !
      !                 clist has to be created by the user program before
      !                 calling this routine. ppm_neighlist_clist should be
      !                 used for this. The reason is that this allows the
      !                 user to have multiple cell lists and call this
      !                 routine for the appropriate one. Do not forget to
      !                 call ppm_neighlist_clist_destroy and to DELLOCATE
      !                 Nm after the cell lists are no longer needed.
      !
      !                 dpd needs to be allocated to proper size before
      !                 calling this routine. Also, this routine is not
      !                 resetting dpd to zero before doing the PP
      !                 interactions. This allows contributions from
      !                 different kernels to be accumulated. If needed,
      !                 set it to zero before calling this routine the 
      !                 first time.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_comp_pp_cell.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/11/04 16:28:53  ivos
      !  Optimized loops for the case of empty cells.
      !
      !  Revision 1.5  2004/10/13 16:42:10  davidch
      !  added support for 3d sph kernels
      !
      !  Revision 1.4  2004/07/26 13:37:39  ivos
      !  bugfix: function pointer argument does not need to be declared
      !  EXTERNAL if there is an explicit interface.
      !
      !  Revision 1.3  2004/07/26 11:46:53  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:45:23  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/07/23 12:57:06  ivos
      !  Preliminary checkin due to module clean-up. Not tested yet, but
      !  it compiles.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_si(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_di(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_sci(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_dci(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_su(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_du(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_scu(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_dcu(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_st(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_cell_dt(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_sct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_cell_dct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nsublist,Nm,clist,cutoff2,dpd,info)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_neighlist
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_neighlist_MkNeighIdx
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: xp
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: pdata
      REAL(MK)   , DIMENSION(:,:), INTENT(  OUT) :: dpd
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
      COMPLEX(MK), DIMENSION(:,:), INTENT(  OUT) :: dpd
#endif
      INTEGER                    , INTENT(IN   ) :: Np,lda,nsublist
#if   __KERNEL == __INTERNAL
      INTEGER                    , INTENT(IN   ) :: kernel
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
#endif
#elif __KERNEL == __LOOKUP_TABLE
      REAL(MK)                   , INTENT(IN   ) :: kpar
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:  ), INTENT(IN   ) :: kernel
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: kernel
#endif
#endif
      LOGICAL                    , INTENT(IN   ) :: lsymm
      TYPE(ppm_type_ptr_to_clist), DIMENSION(:), INTENT(IN) :: clist
      INTEGER    , DIMENSION(:,:), INTENT(IN   ) :: Nm
      REAL(MK)                   , INTENT(IN   ) :: cutoff2
      INTEGER                    , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! counters
      INTEGER                                    :: i,idom,ibox,jbox,idx
      INTEGER                                    :: ipart,jpart,ip,jp
      INTEGER                                    :: cbox,iinter,j,k,ispec
      ! coordinate differences
      REAL(MK)                                   :: dx,dy,dz
      REAL(MK)                                   :: factor,factor2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                                   :: summ,summ2
      REAL(MK)                                   :: eta,dm
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)                                :: summ,summ2
      COMPLEX(MK)                                :: eta,dm
#endif
      ! square of inter particle distance
      REAL(MK)                                   :: dij,dij2,dij4,dij5
      ! start and end particle in a box
      INTEGER                                    :: istart,iend,jstart,jend
      ! cell neighbor lists
      INTEGER, DIMENSION(:,:), POINTER           :: inp,jnp
      ! number of interactions for each cell
      INTEGER, SAVE                              :: nnp
      ! cell offsets for box index
      INTEGER                                    :: n1,n2,nz
      REAL(MK)                                   :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
#if   __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              REAL(MK), DIMENSION(:), INTENT(IN) :: kpar
              REAL(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              COMPLEX(MK), DIMENSION(:), INTENT(IN) :: kpar
              COMPLEX(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#endif
#endif
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_comp_pp_cell',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_comp_pp_cell',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_cell',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_cell',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsublist .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_cell',  &
     &            'nsublist must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (cutoff2 .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_cell',  &
     &            'cutoff2 must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_cell',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Build interaction index lists
      !-------------------------------------------------------------------------
      CALL ppm_neighlist_MkNeighIdx(lsymm,inp,jnp,nnp,info)

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS using symmetry
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          DO idom=1,nsublist
              n1  = Nm(1,idom)
              n2  = Nm(1,idom)*Nm(2,idom)
              nz  = Nm(3,idom)
              IF (ppm_dim .EQ. 2) THEN
                  n2 = 0
                  nz = 2
              ENDIF 
              ! loop over all REAL cells (the -2 in the end does this)
              DO k=0,nz-2
                  DO j=0,Nm(2,idom)-2
                      DO i=0,Nm(1,idom)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                               n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                               n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if box is empty
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
                                      IF (lda .EQ. 1) THEN
                                          summ = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                                          DO jpart=(ipart+1),iend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -   &
     &                                            pdata(1,ip))
                                              summ = summ + dm
                                              dpd(1,jp) = dpd(1,jp) - dm
                                          ENDDO
                                          dpd(1,ip) = dpd(1,ip) + summ
                                      ELSEIF (lda .EQ. 2) THEN
                                          summ = 0.0_MK
                                          summ2 = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                                          DO jpart=(ipart+1),iend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -   &
     &                                            pdata(1,ip))
                                              summ = summ + dm
                                              dpd(1,jp) = dpd(1,jp) - dm
                                              dm   = eta*(pdata(2,jp) -   &
     &                                            pdata(2,ip))
                                              summ2 = summ2 + dm
                                              dpd(2,jp) = dpd(2,jp) - dm
                                          ENDDO
                                          dpd(1,ip) = dpd(1,ip) + summ
                                          dpd(2,ip) = dpd(2,ip) + summ2
                                      ELSE
                                          DO jpart=(ipart+1),iend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              DO ispec=1,lda
                                                  dm = eta*(pdata(ispec,jp)- &
     &                                                pdata(ispec,ip))
                                                  dpd(ispec,ip)=dpd(ispec,ip)+dm
                                                  dpd(ispec,jp)=dpd(ispec,jp)-dm
                                              ENDDO
                                          ENDDO
                                      ENDIF
                                  ENDDO
                              !-------------------------------------------------
                              !  For the other boxes check all particles
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = clist(idom)%lhbx(jbox)
                                  jend   = clist(idom)%lhbx(jbox+1)-1
                                  ! skip this iinter if other box is empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      IF (lda .EQ. 1) THEN
                                          summ = 0.0_MK
                                          ! check against all particles 
                                          ! in the other cell
#ifdef __SXF90
!CDIR NODEP
#endif
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) - &
     &                                            pdata(1,ip))
                                              summ = summ +dm
                                              dpd(1,jp) = dpd(1,jp) - dm
                                          ENDDO
                                          dpd(1,ip) = dpd(1,ip) + summ
                                      ELSEIF (lda .EQ. 2) THEN
                                          summ = 0.0_MK
                                          summ2 = 0.0_MK
                                          ! check against all particles 
                                          ! in the other cell
#ifdef __SXF90
!CDIR NODEP
#endif
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) - &
     &                                            pdata(1,ip))
                                              summ = summ + dm
                                              dpd(1,jp) = dpd(1,jp) - dm
                                              dm   = eta*(pdata(2,jp) - &
     &                                            pdata(2,ip))
                                              summ2 = summ2 + dm
                                              dpd(2,jp) = dpd(2,jp) - dm
                                          ENDDO
                                          dpd(1,ip) = dpd(1,ip) + summ
                                          dpd(2,ip) = dpd(2,ip) + summ2
                                      ELSE
                                          ! check against all particles 
                                          ! in the other cell
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp here... and
                                              ! vice versa to use symmetry.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              DO ispec=1,lda
                                                  dm = eta*(pdata(ispec,jp)- &
     &                                                pdata(ispec,ip))
                                                  dpd(ispec,ip)=dpd(ispec,ip)+dm
                                                  dpd(ispec,jp)=dpd(ispec,jp)-dm
                                              ENDDO
                                          ENDDO
                                      ENDIF
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS not using symmetry
      !-------------------------------------------------------------------------
      ELSE
          DO idom=1,nsublist
              n1  = Nm(1,idom)
              n2  = Nm(1,idom)*Nm(2,idom)
              nz  = Nm(3,idom)
              IF (ppm_dim .EQ. 2) THEN
                  n2 = 0
                  nz = 2
              ENDIF 
              ! loop over all REAL cells (the -2 in the end does this)
              DO k=1,nz-2
                  DO j=1,Nm(2,idom)-2
                      DO i=1,Nm(1,idom)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                               n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                               n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if box is empty
                              !-------------------------------------------------
                              istart = clist(idom)%lhbx(ibox)
                              iend   = clist(idom)%lhbx(ibox+1)-1
                              IF (iend .LT. istart) CYCLE
                              !-------------------------------------------------
                              !  Do all interactions within the box itself
                              !-------------------------------------------------
                              IF (ibox .EQ. jbox) THEN
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      IF (lda .EQ. 1) THEN
                                          DO jpart=istart,iend
                                              jp = clist(idom)%lpdx(jpart)
                                              ! No particle interacts with
                                              ! itself
                                              IF (ip .EQ. jp) CYCLE
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -  &
     &                                            pdata(1,ip))
                                              dpd(1,ip) = dpd(1,ip) + dm
                                          ENDDO
                                      ELSEIF (lda .EQ. 2) THEN
                                          DO jpart=istart,iend
                                              jp = clist(idom)%lpdx(jpart)
                                              ! No particle interacts with
                                              ! itself
                                              IF (ip .EQ. jp) CYCLE
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -  &
     &                                            pdata(1,ip))
                                              dpd(1,ip) = dpd(1,ip) + dm
                                              dm   = eta*(pdata(2,jp) -  &
     &                                            pdata(2,ip))
                                              dpd(2,ip) = dpd(2,ip) + dm
                                          ENDDO
                                      ELSE
                                          DO jpart=istart,iend
                                              jp = clist(idom)%lpdx(jpart)
                                              ! No particle interacts with
                                              ! itself
                                              IF (ip .EQ. jp) CYCLE
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              DO ispec=1,lda
                                                  dm = eta*(pdata(ispec,jp)- &
     &                                                pdata(ispec,ip))
                                                  dpd(ispec,ip)=dpd(ispec,ip)+dm
                                              ENDDO
                                          ENDDO
                                      ENDIF
                                  ENDDO
                              !-------------------------------------------------
                              !  Do interactions with all neighboring boxes
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = clist(idom)%lhbx(jbox)
                                  jend   = clist(idom)%lhbx(jbox+1)-1
                                  ! skip if empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      IF (lda .EQ. 1) THEN
                                          ! check against all particles 
                                          ! in the other cell
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -  &
     &                                            pdata(1,ip))
                                              dpd(1,ip) = dpd(1,ip) + dm
                                          ENDDO
                                      ELSEIF (lda .EQ. 2) THEN
                                          ! check against all particles 
                                          ! in the other cell
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              dm   = eta*(pdata(1,jp) -  &
     &                                            pdata(1,ip))
                                              dpd(1,ip) = dpd(1,ip) + dm
                                              dm   = eta*(pdata(2,jp) -  &
     &                                            pdata(2,ip))
                                              dpd(2,ip) = dpd(2,ip) + dm
                                          ENDDO
                                      ELSE
                                          ! check against all particles 
                                          ! in the other cell
                                          DO jpart=jstart,jend
                                              jp = clist(idom)%lpdx(jpart)
                                              dx  = xp(1,ip) - xp(1,jp)
                                              dy  = xp(2,ip) - xp(2,jp)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz  = xp(3,ip) - xp(3,jp)
                                                  dij = (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij = (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cutoff2) CYCLE
                                              !---------------------------------
                                              ! Particle ip interacts with
                                              ! particle jp.
                                              !---------------------------------
#include "ppm_comp_pp_kernels.inc"
                                              DO ispec=1,lda
                                                  dm = eta*(pdata(ispec,jp)-  &
     &                                                pdata(ispec,ip))
                                                  dpd(ispec,ip)=dpd(ispec,ip)+dm
                                              ENDDO
                                          ENDDO
                                      ENDIF
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_cell',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_dci
#endif
      
#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_cell_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_cell_dct
#endif
#endif
