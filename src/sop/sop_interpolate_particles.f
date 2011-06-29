!!!----------------------------------------------------------------------------!
!!! Interpolate the field variable from the old particles' positions 
!!!      (stored in xp_old) to the new ones (stored in xp)
!!! Assumes that the operators have been computed by&
!!!     correct_diff_operators_interp
!!! 
!!!
!!! Input:
!!!       wp_old:  value at old positions
!!!
!!! Output (all optional, all assummed to be already allocated):
!!!       wp: interpolated value at new positions
!!!       wp_grad: derivative computed with the PSE kernel, using data from old
!!!       positions
!!!       wp_lap: laplacian computed with the PSE kernel, using data from old
!!!       positions
!!!----------------------------------------------------------------------------!

SUBROUTINE sop_interpolate_particles(&
        xp_old,wp_old,xp,Npart,&
        nvlist_cross,vlist_cross,eta,eta_interp,&
        opts,info,wp,wp_grad,wp_lap)

    USE ppm_module_sop_typedef

    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    REAL(MK),     DIMENSION(:,:),           INTENT(IN   )   :: xp_old
    REAL(MK),     DIMENSION(:),             INTENT(IN   )   :: wp_old
    REAL(MK),     DIMENSION(:,:),POINTER,   INTENT(IN   )   :: xp
    INTEGER,                                INTENT(IN   )   :: Npart
    INTEGER,        DIMENSION(:),  POINTER, INTENT(IN   )   :: nvlist_cross
    INTEGER,        DIMENSION(:,:),POINTER, INTENT(IN   )   :: vlist_cross
    REAL(MK),  DIMENSION(:,:),ALLOCATABLE,  INTENT(IN   )   :: eta
    REAL(MK),  DIMENSION(:,:),ALLOCATABLE,  INTENT(IN   )   :: eta_interp
    TYPE(sop_t_opts),  POINTER,             INTENT(IN   )   :: opts
    !!! options
    INTEGER,                                INTENT(  OUT)   :: info

    ! optional arguments
    REAL(MK),DIMENSION(:),POINTER,OPTIONAL,INTENT(INOUT)  :: wp
    REAL(MK),DIMENSION(:,:),POINTER,OPTIONAL,INTENT(INOUT):: wp_grad
    REAL(MK),DIMENSION(:),POINTER,OPTIONAL,INTENT(INOUT)  :: wp_lap

    ! local variables
    INTEGER                         :: ip,iq,ineigh, count
    REAL(KIND(1.d0))                :: t0
    CHARACTER(LEN=256)              :: filename,cbuf
    CHARACTER(LEN=256)              :: caller = 'sop_interpolate_particles'

    !!-------------------------------------------------------------------------!
    !! Initialise
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
    IF (MINVAL(nvlist_cross(1:Npart)) .LT. opts%nneigh_critical) THEN
        WRITE(cbuf,*) 'Too few cross-neighbours, something is wrong'
        CALL pwrite(ppm_rank,caller,cbuf,info)
        info= -1
        GOTO 9999
    ENDIF
#endif
    

    !!-------------------------------------------------------------------------!
    !! Initialise the arrays
    !! and interpolate from old to new
    !!-------------------------------------------------------------------------!
    IF (PRESENT(wp)) THEN
#if debug_verbosity > 0
        WRITE(cbuf,*) 'Computing wp'
        CALL pwrite(ppm_rank,caller,cbuf,info)
#endif
        DO ip = 1,Npart 
            wp(ip) = 0._MK
        ENDDO
        DO ip = 1,Npart
            DO ineigh = 1,nvlist_cross(ip)
                iq = vlist_cross(ineigh,ip)
                wp(ip) = wp(ip) + wp_old(iq) * eta_interp(ineigh,ip)
            ENDDO
        ENDDO
    ELSE
        WRITE(cbuf,*) 'wp should not be an optional argument, after all...'
        !wp is needed to compute wp_grad or wp_lap 
        CALL pwrite(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Initialise the arrays
    !! and compute gradient using data at old particle locations
    !!-------------------------------------------------------------------------!
    IF (PRESENT(wp_grad)) THEN
#if debug_verbosity > 0
        WRITE(cbuf,*) 'Computing wp_grad'
        CALL pwrite(ppm_rank,caller,cbuf,info)
#endif
        DO ip = 1,Npart 
            wp_grad(1:ppm_dim,ip) = 0._MK
        ENDDO
        DO ip = 1,Npart
            DO ineigh = 1,nvlist_cross(ip)
                iq = vlist_cross(ineigh,ip)
                wp_grad(1:ppm_dim,ip) = wp_grad(1:ppm_dim,ip) + &
                    (wp(ip)+wp_old(iq)) * 0.5_MK*       & 
                    (xp_old(1:ppm_dim,iq)-xp(1:ppm_dim,ip)) * eta(ineigh,ip)
            ENDDO
        ENDDO
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Initialise the arrays
    !! and compute laplacian using data at old particle locations
    !!-------------------------------------------------------------------------!
    IF (PRESENT(wp_lap)) THEN
#if debug_verbosity > 0
        WRITE(cbuf,*) 'Computing wp_lap'
        CALL pwrite(ppm_rank,caller,cbuf,info)
#endif
        DO ip = 1,Npart 
            wp_lap(ip) = 0._MK
        ENDDO
        DO ip = 1,Npart
            DO ineigh = 1,nvlist_cross(ip)
                iq = vlist_cross(ineigh,ip)
                wp_lap(ip) = wp_lap(ip) + (wp_old(iq)-wp(ip)) * eta(ineigh,ip) 
            ENDDO
        ENDDO
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999  CONTINUE ! jump here upon error

END SUBROUTINE sop_interpolate_particles


#undef __KIND
