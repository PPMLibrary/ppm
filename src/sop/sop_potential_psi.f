!!-----------------------------------------------------------------------------!
!! Compute interaction potential
!!-----------------------------------------------------------------------------!

!SUBROUTINE sop_potential_psi(xp,D,nvlist,vlist,Npart,Mpart,&
        !Psi_global,Psi_max,info)

SUBROUTINE sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank,ppm_comm,ppm_mpi_kind
    USE ppm_module_io_vtk

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
    INTEGER,                             INTENT(  OUT)   :: info
    REAL(MK),                            INTENT(  OUT)   :: Psi_global
    REAL(MK),                            INTENT(  OUT)   :: Psi_max
    TYPE(sop_t_opts), POINTER,           INTENT(IN   )   :: opts

    ! local variables
    INTEGER                               :: ip,iq,ineigh,iunit,di
    REAL(MK)                              :: dist2,rr,meanD,rd,rc
    REAL(MK),DIMENSION(ppm_dim)           :: dist
    REAL(KIND(1.D0))                      :: t0
    CHARACTER (LEN=256)                   :: caller='sop_potential_psi'
    CHARACTER (LEN=256)                   :: filename,cbuf
    REAL(MK)                              :: Psi_part,attractive_radius
    REAL(MK),DIMENSION(:,:),POINTER       :: xp => NULL()
    REAL(MK),DIMENSION(:  ),POINTER       :: D => NULL()
    REAL(MK),DIMENSION(:  ),POINTER       :: rcp => NULL()
    INTEGER, DIMENSION(:  ),POINTER       :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER       :: vlist => NULL()
    REAL(MK)                              :: rho,coeff
    LOGICAL                               :: no_fusion
    INTEGER,DIMENSION(:),POINTER          :: fuse

#ifdef finnis_sinclair
    REAL(MK)                              :: rho_part
    REAL(MK),PARAMETER                    :: param2_f_s= 1._MK
#endif

    !!-------------------------------------------------------------------------!
    ! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    Psi_global = 0._MK
    !Psi_max = 0._MK

    !!-------------------------------------------------------------------------!
    !! Compute interaction potential
    !!-------------------------------------------------------------------------!

    xp => Get_xp(Particles,with_ghosts=.TRUE.)
    D  => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
    rcp  => Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
    IF (.NOT.Particles%neighlists) THEN
        CALL ppm_write(ppm_rank,caller,&
            'need to compute neighbour lists first',info)
        info = -1
        GOTO 9999
    ENDIF
    nvlist => Particles%nvlist
    vlist => Particles%vlist

    fuse  => Get_wpi(Particles,fuse_id,with_ghosts=.TRUE.)

    particle_loop: DO ip = 1,Particles%Npart
        Psi_part = 0._MK

            !IF (nvlist(ip).GT.30) THEN
        attractive_radius = opts%attractive_radius0
            !ELSE
                !attractive_radius = 0._MK
            !ENDIF

         
        neighbour_loop: DO ineigh = 1,nvlist(ip)
            iq = vlist(ineigh,ip)

            dist = xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq)
            dist2 = SUM(dist**2)
            rr = SQRT(dist2)

#if debug_verbosity > 0
            IF (rr .LE. 1e-12) THEN
                WRITE(cbuf,*) 'Distance between particles too small', &
                    rr,ip,iq,D(ip), D(iq),xp(1:ppm_dim,ip),xp(1:ppm_dim,iq)
                CALL ppm_write(ppm_rank,caller,cbuf,info)
#if debug_verbosity > 1
                WRITE(cbuf,'(A)') 'Part_core_dump'
                CALL ppm_vtk_particle_cloud(cbuf,Particles,info)
                WRITE(cbuf,*) 'Data written in Part_core_dump.vtk'
                CALL ppm_write(ppm_rank,caller,cbuf,info)
#endif
                info = -1
                GOTO 9999
            ENDIF
#endif

            !meanD = 0.5_MK*(D(ip) + D(iq))
            meanD = MIN(D(ip) , D(iq))
            !meanD = MAX(D(ip) , D(iq))

            ! Do not interact with particles which are too far away
            IF (meanD .GT. opts%rcp_over_D * opts%maximum_D) CYCLE

            rd = rr / meanD

            if (fuse(ip)*fuse(iq).GE.1 .and. max(fuse(ip),fuse(iq)).ge.4 ) then 
                no_fusion = .false.
            else
                no_fusion = .true.
            endif

            if (fuse(ip)+fuse(iq).GE.1) then 
                coeff = 1._mk / REAL(MAX(fuse(ip),fuse(iq)),MK)
            else 
                coeff = 1._mk
            endif
            no_fusion = .false.

#include "potential/potential.f90"

        ENDDO neighbour_loop

        Psi_global = Psi_global + Psi_part

    ENDDO particle_loop

    fuse  => set_wpi(Particles,fuse_id,read_only=.TRUE.)

    xp => Set_xp(Particles,read_only=.TRUE.)
    D  => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
    rcp  => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
    nvlist => NULL()
    vlist => NULL()

    !!-------------------------------------------------------------------------!
    !! Compute the global potential of the particles (NOT normalized)
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
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif
    9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_potential_psi

SUBROUTINE sop_plot_potential(opts,filename,info)
    ! write tabulated values of the interaction potential into a file
    ! using parameters from the opts argument

    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    TYPE(sop_t_opts), POINTER,           INTENT(IN   )   :: opts
    CHARACTER(LEN=*),                    INTENT(IN   )   :: filename
    INTEGER,                             INTENT(OUT  )   :: info

    REAL(MK)   :: rd,rho,meanD,coeff
    REAL(MK)   :: Psi_part,attractive_radius
    INTEGER    :: i
    LOGICAL    :: no_fusion


    info = 0
    no_fusion = .FALSE.
    attractive_radius = opts%attractive_radius0
    meanD = 1._mk
    coeff = 1._mk

    OPEN(UNIT=271,FILE=TRIM(ADJUSTL(filename)),IOSTAT=info)
    DO i=1,1000

        Psi_part = 0._mk
        rd = REAL(i,MK)/100._MK

#include "potential/potential.f90"

        WRITE(271,'(2(E22.10,2X))') rd, Psi_part
    ENDDO

    CLOSE(271)

END SUBROUTINE sop_plot_potential

#undef __KIND
