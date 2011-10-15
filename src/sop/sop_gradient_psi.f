!!!----------------------------------------------------------------------------!
!!! Compute gradient of the potential w.r.t particles' positions
!!! Get ghost values for Gradient_Psi
!!!----------------------------------------------------------------------------!
SUBROUTINE sop_gradient_psi(Particles,topo_id,&
        Gradient_Psi,Psi_global,Psi_max,opts,info,gradD)

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank,ppm_comm,ppm_mpi_kind
    USE ppm_module_particles
    USE ppm_module_map_part

    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles),POINTER,       INTENT(INOUT)   :: Particles
    INTEGER,                             INTENT(IN   )   :: topo_id
    REAL(MK),DIMENSION(:,:),POINTER,     INTENT(INOUT)   :: Gradient_Psi
    REAL(MK),                            INTENT(  OUT)   :: Psi_global
    REAL(MK),                            INTENT(  OUT)   :: Psi_max
    TYPE(sop_t_opts), POINTER,           INTENT(IN   )   :: opts
    INTEGER,                             INTENT(  OUT)   :: info

    ! Optional arguments
    REAL(MK),DIMENSION(:,:),POINTER,OPTIONAL,INTENT(IN):: gradD

    ! local variables
    INTEGER                               :: ip,iq,ineigh,iunit,di
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: rr,rr_iso,meanD,nn,rd,rc,scaling_ip,scaling_iq
    REAL(MK),DIMENSION(ppm_dim)           :: dist, grad_vec
    REAL(MK)                              :: Psi_part,gradPsi,attractive_radius
    CHARACTER (LEN = ppm_char)            :: caller='sop_gradient_psi'
    CHARACTER (LEN = ppm_char)            :: filename,cbuf
    REAL(MK),DIMENSION(:,:),POINTER       :: xp => NULL()
    REAL(MK),DIMENSION(:,:),POINTER       :: D => NULL()
    REAL(MK),DIMENSION(:,:),POINTER       :: inv => NULL()
    INTEGER, DIMENSION(:  ),POINTER       :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER       :: vlist => NULL()


    !For debugging
#if debug_verbosity > 1
    REAL(MK),DIMENSION(Particles%Npart)             :: Potential
#endif
    !end debugging

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    DO ip=1,Particles%Mpart
        Gradient_Psi(1:ppm_dim,ip) = 0._MK
    ENDDO

    Psi_global = 0._MK
    Psi_max   = 0._MK

    Particles%h_min     = HUGE(1._MK)

    xp => Get_xp(Particles,with_ghosts=.TRUE.)
    D  => Get_wpv(Particles,Particles%D_id,with_ghosts=.TRUE.)
    inv  => Get_wpv(Particles,Particles%G_id,with_ghosts=.TRUE.)
    IF (.NOT.Particles%neighlists) THEN
        CALL ppm_write(ppm_rank,caller,&
            'need to compute neighbour lists first',info)
        info = -1
        GOTO 9999
    ENDIF
    nvlist => Particles%nvlist
    vlist => Particles%vlist

   ! HAECKIC: combine aniso with iso

    !!-------------------------------------------------------------------------!
    !! Compute interaction potential and its gradient
    !!-------------------------------------------------------------------------!
    particle_loop: DO ip = 1,Particles%Npart
        Psi_part = 0._MK
        
        !
        nn = HUGE(1._MK)
#ifdef finnis_sinclair
        rho_part = 0._MK
#endif

        !Particles with few neighbours are a bit more reluctant to fuse
        !(not really necessary, but makes insertion/deletion a bit faster
        !in some cases)
        !attractive_radius = attractive_radius0 * &
         !MIN(REAL(nvlist(ip)-22,MK)/10._MK,1._MK)

         ! haeckic: keep that??
!          IF (nvlist(ip).GT.30) THEN
                attractive_radius = opts%attractive_radius0
!          ELSE
!                attractive_radius = 0._MK
!          ENDIF

        neighbour_loop: DO ineigh = 1,nvlist(ip)
            iq = vlist(ineigh,ip)

            ! rr: distance in unit circle
            ! dist: vector from current point to neighbor point: xq - xp
            ! rr_iso: distance in real space

            CALL particles_anisotropic_distance(Particles,ip,iq,Particles%G_id,rr,info)
            dist = xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq)

            IF (rr .LE. 1e-12) THEN
                WRITE(cbuf,*) 'Distance between particles less than 1e-12! ',&
                rr, ' ',ip,' ',iq,' xp = ',xp(1,ip),xp(2,ip),&
                    ' xq = ',xp(1,iq),xp(2,iq)
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                info = -1
                GOTO 9999
            ENDIF

            IF (ppm_dim .EQ. 2) THEN
               grad_vec(1) = - inv(1,ip)*(inv(1,ip)*dist(1) + inv(2,ip)*dist(2)) &
               &             - inv(3,ip)*(inv(3,ip)*dist(1) + inv(4,ip)*dist(2)) 
               grad_vec(2) = - inv(2,ip)*(inv(1,ip)*dist(1) + inv(2,ip)*dist(2)) &
               &             - inv(4,ip)*(inv(3,ip)*dist(1) + inv(4,ip)*dist(2))
            ELSE
               grad_vec(1) = - inv(1,ip)*(inv(1,ip)*dist(1) + inv(2,ip)*dist(2) + inv(3,ip)*dist(3)) &
               &             - inv(4,ip)*(inv(4,ip)*dist(1) + inv(5,ip)*dist(2) + inv(6,ip)*dist(3)) &
               &             - inv(7,ip)*(inv(7,ip)*dist(1) + inv(8,ip)*dist(2) + inv(9,ip)*dist(3)) 
               grad_vec(2) = - inv(2,ip)*(inv(1,ip)*dist(1) + inv(2,ip)*dist(2) + inv(3,ip)*dist(3)) &
               &             - inv(5,ip)*(inv(4,ip)*dist(1) + inv(5,ip)*dist(2) + inv(6,ip)*dist(3)) &
               &             - inv(8,ip)*(inv(7,ip)*dist(1) + inv(8,ip)*dist(2) + inv(9,ip)*dist(3))
               grad_vec(3) = - inv(3,ip)*(inv(1,ip)*dist(1) + inv(2,ip)*dist(2) + inv(3,ip)*dist(3)) &
               &             - inv(6,ip)*(inv(4,ip)*dist(1) + inv(5,ip)*dist(2) + inv(6,ip)*dist(3)) &
               &             - inv(9,ip)*(inv(7,ip)*dist(1) + inv(8,ip)*dist(2) + inv(9,ip)*dist(3))
            ENDIF

            CALL particles_shorter_axis(Particles,ip,Particles%G_id,scaling_ip,info)
            scaling_ip = scaling_ip**2

            
            ! HAECKIC: The variables
            ! rr: distance used in gradient potential
            ! scaling_ip: scaling used for ip
            ! scaling_iq: scaling used for iq
            ! grad_vec: additional term for gradient in anisotropic case
            
            !meanD = 0.5_MK*(D(ip) + D(iq))
            !meanD = MAX(D(ip) , D(iq))
            !meanD = MIN(D(ip) , D(iq))

            ! Do not interact with particles which are too far away
            ! haeckic: do we need that?
            !IF (meanD .GT. opts%rcp_over_D * opts%maximum_D) CYCLE

            !------------------------------------------------------------------!
            !Compute the gradients with respect to rpq
            !------------------------------------------------------------------!
            !------------------------------------------------------------------!
            ! here we can choose between different interaction potentials
            !------------------------------------------------------------------!
!             rd = rr / meanD
!             rc = rr / (MIN(rcp(ip),rcp(iq)))

!haeckic: check the way the potential is used

#include "potential/potential_gradient.f90"



            ! add part to ip and subtract part from iq
!             write(*,*) 'Particle', ip, scaling_ip, rr, grad_vec, gradPsi
!             write(*,*) 'positions', xp(1:ppm_dim,ip), xp(1:ppm_dim,iq), rr_iso
            Gradient_Psi(1:ppm_dim,ip)=Gradient_Psi(1:ppm_dim,ip) &
                  + scaling_ip/rr * grad_vec * gradPsi
            Gradient_Psi(1:ppm_dim,iq)=Gradient_Psi(1:ppm_dim,iq) &
                  - scaling_ip/rr * grad_vec * gradPsi

        ENDDO neighbour_loop

        Psi_global = Psi_global + Psi_part

        !For debugging, not needed
#if debug_verbosity > 1
        Potential(ip)=Psi_part
#endif
        !end debugging

        !!---------------------------------------------------------------------!
        !! Cap the values of the gradient such that particle displacements 
        !! cannot be larger than the size of the ghost layers.
        !! FIXME (account for variable ghost layer sizes)
        !! FIXME vectorize loops
        !!---------------------------------------------------------------------!
        ! HAECKIC: the capping is done after the particles loop, as it is asymmetric

        !----------------------------------------------------------------------!
        ! Save the minimum and maximum separation distance on this processor
        ! (used for checking stability conditions, for example)
        !----------------------------------------------------------------------!
        ! haeckic: keep that?
        Particles%h_min = MIN(Particles%h_min,nn)

#if debug_verbosity > 1
        IF (ISNAN(Psi_part)) THEN
            WRITE(*,*) 'Psi_part = ',Psi_part, ip, gradPsi, D(ip), nn
            WRITE(*,*) xp(1:ppm_dim,ip)
            WRITE(*,*) 'Grad = ',Gradient_Psi(1:ppm_dim,ip)
            DO ineigh=1,nvlist(ip)
                iq = vlist(ineigh,ip)
                WRITE(*,*) 'xp(iq)= ',xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq)
            ENDDO
            DO ineigh=1,nvlist(ip)
                iq = vlist(ineigh,ip)
                WRITE(*,*) 'D(iq)= ',D(iq)
            ENDDO
            info = -1
            GOTO 9999
        ENDIF
#endif

    ENDDO particle_loop

   ! HAECKIC: scaling is done after loop
   DO ip = 1,Particles%Npart
         ! The scaling needs to be done after the particle loop, because asymmetric
         ! scaled to real distance, not the transformed one
         ! If the gradient length is bigger than distance to nearest neighbor?
         CALL particles_shorter_axis(Particles,ip,Particles%G_id,nn,info)
         IF (SUM(Gradient_Psi(1:ppm_dim,ip)**2) .GT. MIN(nn**2,Particles%cutoff**2)) THEN
               Gradient_Psi(1:ppm_dim,ip) = Gradient_Psi(1:ppm_dim,ip) * &
                  MIN(nn,Particles%cutoff) * 0.9_MK / SQRT(SUM(Gradient_Psi(1:ppm_dim,ip)**2))
         ENDIF
   ENDDO 

    xp => Set_xp(Particles,read_only=.TRUE.)
    D  => Set_wpv(Particles,Particles%D_id,read_only=.TRUE.)
    inv  => Set_wpv(Particles,Particles%G_id,read_only=.TRUE.)
    nvlist => NULL()
    vlist => NULL()


#if debug_verbosity > 2
    !For debugging, not needed
    CALL sop_dump_debug(Potential,Particles%Npart,477,info)
    !end debugging
#endif

    !!-------------------------------------------------------------------------!
    !! Compute the global potential of the particles 
    !! (it can be viewed as an objective function)
    !!-------------------------------------------------------------------------!

#ifdef __MPI
    CALL MPI_Allreduce(Psi_global,Psi_global,1,ppm_mpi_kind,MPI_SUM,ppm_comm,info)
    CALL MPI_Allreduce(Psi_max,Psi_max,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'MPI_Allreduce failed',info)
        info = -1
        GOTO 9999
    ENDIF
#endif

    !!-------------------------------------------------------------------------!
    !! Get ghosts for gradient_psi (then, during the linesearch, 
    !! we can move ghost particles directly, without communicating)
    !!-------------------------------------------------------------------------!

    CALL ppm_map_part_push(Gradient_Psi,ppm_dim,Particles%Npart,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'ppm_map_part_push failed',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL ppm_map_part_send(Particles%Npart,Particles%Mpart,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'ppm_map_part_send failed',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL ppm_map_part_pop(Gradient_Psi,ppm_dim,Particles%Npart, &
        Particles%Mpart,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'ppm_map_part_pop failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_gradient_psi


#undef __KIND
