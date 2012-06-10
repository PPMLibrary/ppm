!!!----------------------------------------------------------------------------!
!!! Approximate values of wp at new particle positions from the values wp_old
!!! at xp_old by doing a weighted average (first order interpolation)
!!!
!!! and get the ghosts
!!!
!!! Not very accurate, but much faster than computing corrected interpolation
!!! kernels.
!!!
!!! Input:
!!!       wp_old:  value at old positions
!!!
!!! Output:  
!!!       wp: approximate values at new positions, inc. ghost particles
!!!
!!! Remark: - wp is assumed to be already allocated
!!!         - Requires local communications (get_ghosts) 
!!!----------------------------------------------------------------------------!

#if   __LDA == 1
SUBROUTINE DTYPE(sop_approx_wp_1d)(xp_old,wp_old,scale_old,xp,wp,&
        Npart,Mpart,nvlist_cross,vlist_cross,info)
#elif __LDA == 2
SUBROUTINE DTYPE(sop_approx_wp_2d)(xp_old,wp_old,lda,scale_old,xp,wp,&
        Npart,Mpart,nvlist_cross,vlist_cross,info)
#endif

    USE ppm_module_map_part
    IMPLICIT NONE

#ifdef __MPI
    INCLUDE 'mpif.h'
#endif

    DEFINE_MK()
    ! arguments
    REAL(MK),     DIMENSION(:,:),         INTENT(IN   )   :: xp_old
    REAL(MK),     DIMENSION(:,:),POINTER, INTENT(IN   )   :: xp
#if   __LDA == 1
    REAL(MK),     DIMENSION(:),           INTENT(IN   )   :: wp_old
    REAL(MK),     DIMENSION(:),  POINTER, INTENT(INOUT)   :: wp
#elif __LDA == 2
    REAL(MK),     DIMENSION(:,:),         INTENT(IN   )   :: wp_old
    REAL(MK),     DIMENSION(:,:),POINTER, INTENT(INOUT)   :: wp
    INTEGER,                              INTENT(IN   )   :: lda
#endif
    REAL(MK),     DIMENSION(:),           INTENT(IN   )   :: scale_old
    INTEGER,                              INTENT(IN   )   :: Npart
    INTEGER,                              INTENT(INOUT)   :: Mpart
    INTEGER,      DIMENSION(:),  POINTER, INTENT(IN   )   :: nvlist_cross
    INTEGER,      DIMENSION(:,:),POINTER, INTENT(IN   )   :: vlist_cross
    INTEGER,                              INTENT(  OUT)   :: info


        ! local variables
        INTEGER                         :: ip,iq,ineigh, count
        REAL(KIND(1.d0))                :: t0
        CHARACTER(LEN=256)              :: filename,cbuf
#if   __LDA == 1
        CHARACTER(LEN=256)              :: caller = 'sop_approx_wp_1d'
#elif __LDA == 2
        CHARACTER(LEN=256)              :: caller = 'sop_approx_wp_2d'
#endif
        REAL(MK)                      :: weight,weight_sum
        REAL(MK)                      :: almostzero,h_local

        !!---------------------------------------------------------------------!
        !! Initialise
        !!---------------------------------------------------------------------!
        info = 0
#if debug_verbosity > 0
        CALL substart(caller,t0,info)
#endif
        almostzero = 100._MK*EPSILON(1._MK)

        !!---------------------------------------------------------------------!
        !! Approximate wp as the weighted sum of wp_old
        !! The weights are taken as the squared distances beween xp and xp_old
        !! EDIT: take the 4th power instead.
        !!---------------------------------------------------------------------!
        DO ip = 1,Npart 
#if   __LDA == 1
            wp(ip) = 0._MK
#elif __LDA == 2
            wp(1:lda,ip) = 0._MK
#endif
        ENDDO
        
        particle_loop: DO ip = 1,Npart

            IF( nvlist_cross(ip) .EQ. 0) THEN
                !CYCLE particle_loop ! D is not updated for this particle
                WRITE(cbuf,*) 'error in sop_approx_D.',&
                    ' Cross lists should not be empty'
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                info = -1
                GOTO 9999
            ENDIF

            weight_sum = 0._MK

            h_local = 0._MK
            DO ineigh = 1,nvlist_cross(ip)
                h_local = h_local + scale_old(vlist_cross(ineigh,ip))
            ENDDO
            h_local = h_local / nvlist_cross(ip)
            neighbour_loop: DO ineigh = 1,nvlist_cross(ip)

                iq = vlist_cross(ineigh,ip)
                weight = SUM(((xp_old(1:ppm_dim,iq)-xp(1:ppm_dim,ip))/h_local)**2) ** 2
                IF (weight .LT. almostzero) THEN
#if   __LDA == 1
                    wp(ip) = wp_old(iq)
#elif __LDA == 2
                    wp(1:lda,ip) = wp_old(1:lda,iq)
#endif
                    CYCLE particle_loop
                ELSE
                    weight_sum = weight_sum + 1._MK / weight
#if   __LDA == 1
                    wp(ip) = wp(ip) + wp_old(iq) / weight
#elif __LDA == 2
                    wp(1:lda,ip) = wp(1:lda,ip) + wp_old(1:lda,iq) / weight
#endif
                ENDIF

            ENDDO neighbour_loop

#if   __LDA == 1
            wp(ip) = wp(ip) / weight_sum
#elif __LDA == 2
            wp(1:lda,ip) = wp(1:lda,ip) / weight_sum
#endif

        ENDDO particle_loop

        !!---------------------------------------------------------------------!
        !! Get ghost values
        !!---------------------------------------------------------------------!
#if   __LDA == 1
        CALL ppm_map_part_push(wp,Npart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (push) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_send(Npart,Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (send) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_pop(wp,Npart,Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (pop) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
#elif __LDA == 2
        CALL ppm_map_part_push(wp,lda,Npart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (push) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_send(Npart,Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (send) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_pop(wp,lda,Npart,Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (pop) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
#endif

        !!---------------------------------------------------------------------!
        !! Finalize
        !!---------------------------------------------------------------------!
#if debug_verbosity > 0
        CALL substop(caller,t0,info)
#endif

        9999  CONTINUE ! jump here upon error

#if   __LDA == 1
    END SUBROUTINE DTYPE(sop_approx_wp_1d)
#elif __LDA == 2
END SUBROUTINE DTYPE(sop_approx_wp_2d)
#endif

#undef __LDA
