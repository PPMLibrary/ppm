      !!-----------------------------------------------------------------------------!
      !! Compute interaction potential
      !!-----------------------------------------------------------------------------!

      !SUBROUTINE sop_potential_psi(xp,D,nvlist,vlist,Npart,Mpart,&
              !Psi_global,Psi_max,info)

      SUBROUTINE DTYPE(sop_potential_psi)(Particles,Psi_global,Psi_max,opts,info)

          USE ppm_module_mpi
          USE ppm_module_data, ONLY: ppm_dim,ppm_rank,ppm_comm,ppm_mpi_kind
          USE ppm_module_io_vtk
          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          TYPE(DTYPE(ppm_t_particles)),POINTER,INTENT(INOUT)   :: Particles
          INTEGER,                             INTENT(  OUT)   :: info
          REAL(MK),                            INTENT(  OUT)   :: Psi_global
          REAL(MK),                            INTENT(  OUT)   :: Psi_max
          TYPE(DTYPE(sop_t_opts)), POINTER,    INTENT(IN   )   :: opts

          ! local variables
          INTEGER                               :: ip,iq,ineigh,iunit,di
          REAL(MK)                              :: rr,meanD,rd,rc
          REAL(ppm_kind_double) :: t0
          CHARACTER (LEN=256)                   :: caller='sop_potential_psi'
          CHARACTER (LEN=256)                   :: filename,cbuf
          REAL(MK)                              :: Psi_part,attractive_radius
          REAL(MK),DIMENSION(:,:),POINTER       :: xp => NULL()
          REAL(MK),DIMENSION(:  ),POINTER       :: D => NULL()
          REAL(MK),DIMENSION(:  ),POINTER       :: rcp => NULL()
          INTEGER, DIMENSION(:  ),POINTER       :: nvlist => NULL()
          INTEGER, DIMENSION(:,:),POINTER       :: vlist => NULL()
          REAL(MK)                              :: rho,coeff,Psi_at_cutoff
          LOGICAL                               :: no_fusion
          INTEGER,DIMENSION(:),POINTER          :: fuse

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

          !offset potential so that it is zero at r=r_cutoff
          rho = opts%param_morse
          Psi_at_cutoff = (-rho**(-4._MK*opts%rcp_over_D) + &
              0.8_MK*rho**(1.0_MK-5._MK*opts%rcp_over_D))
          coeff = 1.0_MK


          attractive_radius = opts%attractive_radius0
          particle_loop: DO ip = 1,Particles%Npart
              Psi_part = 0._MK

              neighbour_loop: DO ineigh = 1,nvlist(ip)
                  iq = vlist(ineigh,ip)

                  rr = SQRT(SUM((xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq))**2))

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

                  meanD = MIN(D(ip),D(iq))

                  rd = rr / meanD

                  !if (fuse(ip)*fuse(iq).GE.1 .and. max(fuse(ip),fuse(iq)).ge.4 ) then
                  if (max(fuse(ip),fuse(iq)).ge.4 ) then
                      no_fusion = .FALSE.
                  else
                      no_fusion = .TRUE.
                  endif

                  !if (fuse(ip)+fuse(iq).GE.1) then
                      !coeff = 1.0_MK / REAL(MAX(fuse(ip),fuse(iq)),MK)
                  !else
                      !coeff = 1.0_MK
                  !endif

                  !------------------------------------------------------------------!
                  ! here we can choose between different interaction potentials
                  !------------------------------------------------------------------!
#if   __SOP_POTENTIAL == __MORSE
#include "potential_morse/potential.f90"
#elif __SOP_POTENTIAL == __SCHRADER
#include "potential_schrader/potential.f90"
#elif __SOP_POTENTIAL == __REPULSIVE
#include "potential_repulsive/potential.f90"
#else
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_sub_failed,caller,&
                  'pair-potential not defined - compile with -D_SOP_POTENTIAL',&
              __LINE__,info)
#endif

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

      END SUBROUTINE DTYPE(sop_potential_psi)

      SUBROUTINE DTYPE(sop_plot_potential)(opts,filename,info)
          ! write tabulated values of the interaction potential into a file
          ! using parameters from the opts argument

          USE ppm_module_mpi
          IMPLICIT NONE

          DEFINE_MK()
          TYPE(DTYPE(sop_t_opts)), POINTER,    INTENT(IN   )   :: opts
          CHARACTER(LEN=*),                    INTENT(IN   )   :: filename
          INTEGER,                             INTENT(OUT  )   :: info

          REAL(MK)   :: rd,rho,meanD,coeff,Psi_at_cutoff
          REAL(MK)   :: Psi_part,attractive_radius,Psi_max
          INTEGER    :: i
          LOGICAL    :: no_fusion


          info = 0
          no_fusion = .FALSE.
          attractive_radius = opts%attractive_radius0
          meanD = 1.0_MK
          coeff = 1.0_MK
          rho = opts%param_morse
          Psi_at_cutoff = (-rho**(-4._MK*opts%rcp_over_D) + &
              0.8_MK*rho**(1.0_MK-5._MK*opts%rcp_over_D))

          OPEN(UNIT=271,FILE=TRIM(ADJUSTL(filename)),IOSTAT=info)
          DO i=1,1000

              Psi_part = 0.0_MK
              rd = REAL(i,MK)/100._MK

#if   __SOP_POTENTIAL == __MORSE
#include "potential_morse/potential.f90"
#elif __SOP_POTENTIAL == __SCHRADER
#include "potential_schrader/potential.f90"
#elif __SOP_POTENTIAL == __REPULSIVE
#include "potential_repulsive/potential.f90"
#else
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_sub_failed,caller,&
                  'pair-potential not defined - compile with -D_SOP_POTENTIAL',&
              __LINE__,info)
#endif

              WRITE(271,'(2(E22.10,2X))') rd, Psi_part
          ENDDO

          CLOSE(271)

      END SUBROUTINE DTYPE(sop_plot_potential)
