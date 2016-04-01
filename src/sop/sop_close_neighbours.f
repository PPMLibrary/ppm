      SUBROUTINE DTYPE(sop_close_neighbours)(Particles,opts,info)
          !!! Compute and store the number of neighbours that each
          !!! particle has in a circle of radius D_tilde

          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          TYPE(DTYPE(ppm_t_particles)), POINTER,INTENT(INOUT)   :: Particles
          TYPE(DTYPE(sop_t_opts)), POINTER,     INTENT(IN   )   :: opts
          INTEGER,                              INTENT(  OUT)   :: info
          !local variables

          REAL(MK),DIMENSION(:,:), POINTER      :: xp
          INTEGER,DIMENSION(:),    POINTER      :: nvlist
          INTEGER,DIMENSION(:,:),  POINTER      :: vlist
          INTEGER                               :: ip,iq,ineigh
          REAL(MK),DIMENSION(:), POINTER        :: D
          REAL(MK),DIMENSION(:), POINTER        :: Dtilde
          REAL(MK)                              :: rr
          REAL(MK)                              :: nn,max_nn,avg_nn
          INTEGER                               :: close_neigh
          INTEGER                               :: very_close_neigh
          CHARACTER(LEN=ppm_char)               :: cbuf
          CHARACTER(LEN=ppm_char)               :: caller = 'sop_close_neighbours'
          INTEGER,DIMENSION(:),POINTER          :: nb_neigh
          INTEGER                               :: nb_close_theo, nb_fuse_neigh

          info = 0
          nvlist => Particles%nvlist
          vlist => Particles%vlist
          xp => Get_xp(Particles,with_ghosts=.TRUE.)
          D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)

          Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)

          nb_neigh => Get_wpi(Particles,nb_neigh_id)


          IF (ppm_dim.EQ.2) THEN
             nb_close_theo = 6
          ELSE
             nb_close_theo = 6
          ENDIF

          particle_loop: DO ip = 1,Particles%Npart
              close_neigh = 0

              neighbour_loop: DO ineigh = 1,nvlist(ip)
                  iq = vlist(ineigh,ip)

                  rr = SQRT(SUM((xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq))**2)) / D(ip)
                  !Dtilde(ip)
                  !MIN(Dtilde(ip),Dtilde(iq))
                  !MIN(D(ip),D(iq))

                  !compute nearest-neighbour distance for each particle
                  ! rescaled by the MIN of D(ip) and D(iq)
                  IF (rr .LT. 1.0_mk) THEN
                      close_neigh = close_neigh + 1

                      !sort vlist (only for the close neighbours, and only
                      ! needed by fuse2_particles
                      vlist(ineigh,ip) = vlist(close_neigh,ip)
                      vlist(close_neigh,ip) = iq
                  ENDIF

              ENDDO neighbour_loop

              IF (close_neigh .LE. nb_close_theo-2) THEN
                 adaptation_ok = .false.
              ENDIF
              !IF (close_neigh .GE. 2*nb_close_theo) THEN
                      !adaptation_ok = .false.
              !ENDIF

              nb_neigh(ip) = close_neigh


          ENDDO particle_loop

          xp => Set_xp(Particles,read_only=.TRUE.)
          D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
          Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)
          nb_neigh => Set_wpi(Particles,nb_neigh_id)

      END SUBROUTINE DTYPE(sop_close_neighbours)

