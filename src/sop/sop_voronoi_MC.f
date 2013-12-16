      !!-----------------------------------------------------------------------------!
      !! Compute voronoi volumes and c.o.m using Monte Carlo sampling
      !!-----------------------------------------------------------------------------!
      SUBROUTINE DTYPE(sop_voronoi_MC)(Particles,opts,info)

          USE ppm_module_data, ONLY: ppm_dim,ppm_rank,ppm_comm,ppm_mpi_kind
          USE ppm_module_io_vtk
          USE ppm_module_kdtree

          IMPLICIT NONE
#ifdef __MPI
          INCLUDE 'mpif.h'
#endif

          DEFINE_MK()
          ! arguments
          TYPE(DTYPE(ppm_t_particles)),POINTER,INTENT(INOUT)   :: Particles
          TYPE(DTYPE(sop_t_opts)), POINTER,    INTENT(IN   )   :: opts
          INTEGER,                             INTENT(  OUT)   :: info

          ! local variables
          INTEGER                               :: ip,iq,ineigh,iunit,di,j
          REAL(MK)                              :: rr,meanD,rd,rc
          REAL(KIND(1.D0))                      :: t0
          CHARACTER (LEN=256)                   :: caller='sop_voronoi_MC'
          CHARACTER (LEN=256)                   :: filename,cbuf
          REAL(MK),DIMENSION(:,:),POINTER       :: xp => NULL()
          REAL(MK),DIMENSION(:  ),POINTER       :: D => NULL()
          REAL(MK),DIMENSION(:  ),POINTER       :: rcp => NULL()
          INTEGER, DIMENSION(:  ),POINTER       :: nvlist => NULL()
          INTEGER, DIMENSION(:,:),POINTER       :: vlist => NULL()

          REAL(MK),DIMENSION(:  ),POINTER       :: vor => NULL()
          INTEGER, PARAMETER                    :: nsamples = 50
          REAL(MK), PARAMETER                   :: twopi = 2._MK*ACOS(-1._MK)
          REAL(MK), PARAMETER                   :: pi = ACOS(-1._MK)
          REAL(MK),DIMENSION(:,:,:),POINTER     :: xsamples
          REAL(MK)                              :: x_s,y_s

          TYPE(ppm_t_topo), POINTER             :: topo
          TYPE(DTYPE(kdtree2)),POINTER            :: tree
          TYPE(DTYPE(kdtree2_result)),ALLOCATABLE :: results(:)
          INTEGER                               :: knn = 1
          !!-------------------------------------------------------------------------!
          ! Initialize
          !!-------------------------------------------------------------------------!
          info = 0
#if debug_verbosity > 0
          CALL substart(caller,t0,info)
#endif

          IF (.NOT.Particles%neighlists) THEN
              CALL ppm_write(ppm_rank,caller,&
                  'need to compute neighbour lists first',info)
              info = -1
              GOTO 9999
          ENDIF


          CALL particles_allocate_wps(Particles,voronoi_id,info,&
              iopt=ppm_param_alloc_fit,name='voronoi volume')
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,caller,&
                  'particles_allocate_wps failed', __LINE__,info)
              GOTO 9999
          ENDIF

          !!-------------------------------------------------------------------------!
          !! Generate MC samples
          !!-------------------------------------------------------------------------!
          ldc(1) = ppm_dim
          ldc(2) = Particles%Npart
          ldc(3) = nsamples
          CALL ppm_alloc(xsamples,ldc,ppm_param_alloc_fit,info)

          xp => Get_xp(Particles,with_ghosts=.TRUE.)
          rcp  => Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
          D  => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
          vor  => Get_wps(Particles,voronoi_id)

          CALL random_number(xsamples)
          DO ip=1,Particles%Npart
              DO j=1,nsamples
                  x_s = sqrt(D(ip)*xsamples(1,ip,j))*cos(twopi*xsamples(2,ip,j))
                  y_s = sqrt(D(ip)*xsamples(1,ip,j))*sin(twopi*xsamples(2,ip,j))
                  xsamples(1,ip,j)= x_s + xp(1,ip)
                  xsamples(2,ip,j)= y_s + xp(2,ip)
              ENDDO
          ENDDO

          tree => kdtree2_create(Particles%xp(1:ppm_dim,1:Particles%Mpart),&
              sort=.true.,rearrange=.true.)
          allocate(results(knn+1))

          !Find nearest particle of each query point
          DO ip=1,Particles%Npart
              vor(ip) = 0._MK
              DO j=1,nsamples
                  call kdtree2_n_nearest(tp=tree,qv=xsamples(1:ppm_dim,ip,j),&
                      nn=knn+1,results=results)

                  IF (results(1)%idx .EQ. ip) THEN
                      vor(ip) = vor(ip) + 1._MK
                  ENDIF
              ENDDO
              vor(ip) = vor(ip)/REAL(nsamples,MK) ! * pi*D(ip)**2
          ENDDO

          call kdtree2_destroy(tree)
          deallocate(results)


          xp => Set_xp(Particles,read_only=.TRUE.)
          D  => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
          rcp  => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
          vor  => Set_wps(Particles,voronoi_id)
          nvlist => NULL()
          vlist => NULL()


          !!-------------------------------------------------------------------------!
          !! Finalize
          !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
          CALL substop(caller,t0,info)
#endif
          9999 CONTINUE ! jump here upon error


      END SUBROUTINE DTYPE(sop_voronoi_MC)
