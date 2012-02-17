         SUBROUTINE DTYPE(ppm_vtk_particle_cloud)(filename, Particles, info, &
              step, with_ghosts, with_nvlist, wpi_list, wps_list,     &
              wpv_list, wpv_field_list)
           !--------------------------------------------------------------------
           !  Arguments
           !--------------------------------------------------------------------
           DEFINE_MK()
           CHARACTER(LEN=*),                INTENT(IN   ) :: filename
           TYPE(DTYPE(ppm_t_particles)),POINTER,INTENT(IN):: Particles
           INTEGER,                         INTENT(  OUT) :: info
           INTEGER,               OPTIONAL, INTENT(IN   ) :: step
           LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_ghosts
           LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_nvlist
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wpi_list
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wps_list
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wpv_list
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wpv_field_list
           !--------------------------------------------------------------------
           !  Variables
           !--------------------------------------------------------------------
           CHARACTER(LEN=*), PARAMETER :: caller='ppm_vtk_particle_cloud'
           CHARACTER(LEN=ppm_char)              :: scratch
           CHARACTER(LEN=ppm_char)              :: fname
           REAL(ppm_kind_double)                :: t0
           INTEGER                              :: i, j, k, l, nd, N
           INTEGER                              :: nb_wps, nb_wpv, nb_wpv_field
           INTEGER                              :: nb_wpi
           INTEGER, DIMENSION(:),   ALLOCATABLE :: wpi_l, wps_l, wpv_l, wpv_field_l
           LOGICAL                              :: ghosts
           LOGICAL                              :: nvlist
           REAL(MK), DIMENSION(:,:),POINTER     :: xp  => NULL()
           REAL(MK), DIMENSION(:),  POINTER     :: wp  => NULL()
           INTEGER, DIMENSION(:),   POINTER     :: wpi  => NULL()
           !--------------------------------------------------------------------
           !  Code
           !--------------------------------------------------------------------
           CALL substart(caller,t0,info)
           
           IF (PRESENT(step)) THEN
              WRITE(fname,'(A,A,I0)') &
                   filename(1:LEN_TRIM(filename)), '.', step
           ELSE
              fname = filename
           END IF

           ! print ghosts?
           IF (PRESENT(with_ghosts)) THEN
              ghosts = with_ghosts
           ELSE
              ghosts = .FALSE.
           END IF

           ! print nvlists?
           IF (PRESENT(with_nvlist)) THEN
              nvlist = with_nvlist
           ELSE
              nvlist = .FALSE.
           END IF


           IF (nvlist .AND. ghosts) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,caller,   &
                   &  'printout of nvlist for ghosts not supported (yet)',&
                   &  __LINE__,info)
               GOTO 9999
           ENDIF

           ! create the list of properties to print
           ! integer property
           IF(PRESENT(wpi_list)) THEN
              nb_wpi=SIZE(wpi_list)
              ALLOCATE(wpi_l(nb_wpi),STAT=info)
              wpi_l=wpi_list
              DO i=1,nb_wpi
                 IF (wpi_l(i).GT.Particles%max_wpiid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'integer property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (.NOT.Particles%wpi(wpi_l(i))%is_mapped) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              !printout all properties i that are mapped
              nb_wpi = 0
              DO i=1,Particles%max_wpiid
                 IF (Particles%wpi(i)%is_mapped) &
                      nb_wpi = nb_wpi + 1
              ENDDO
              ALLOCATE(wpi_l(nb_wpi),STAT=info)
              nb_wpi = 0
              DO i=1,Particles%max_wpiid
                 IF (Particles%wpi(i)%is_mapped) THEN
                    nb_wpi = nb_wpi + 1
                    wpi_l(nb_wpi) = i
                 ENDIF
              ENDDO
           ENDIF

           ! scalar property
           IF(PRESENT(wps_list)) THEN
              nb_wps=SIZE(wps_list)
              ALLOCATE(wps_l(nb_wps),STAT=info)
              wps_l=wps_list
              DO i=1,nb_wps
                 IF (wps_l(i).GT.Particles%max_wpsid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'scalar property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (.NOT.Particles%wps(wps_l(i))%is_mapped) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              !printout all properties i that are mapped
              nb_wps = 0
              DO i=1,Particles%max_wpsid
                 IF (Particles%wps(i)%is_mapped) &
                      nb_wps = nb_wps + 1
              ENDDO
              ALLOCATE(wps_l(nb_wps),STAT=info)
              nb_wps = 0
              DO i=1,Particles%max_wpsid
                 IF (Particles%wps(i)%is_mapped) THEN
                    nb_wps = nb_wps + 1
                    wps_l(nb_wps) = i
                 ENDIF
              ENDDO
           ENDIF

           ! vector property
           IF(PRESENT(wpv_list)) THEN
              nb_wpv=SIZE(wpv_list)
              ALLOCATE(wpv_l(nb_wpv),STAT=info)
              wpv_l=wpv_list
              DO i=1,nb_wpv
                 IF (wpv_l(i).GT.Particles%max_wpvid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'vector property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (.NOT.Particles%wpv(wpv_l(i))%is_mapped) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              !printout all properties i that are mapped
              nb_wpv = 0
              DO i=1,Particles%max_wpvid
                 IF (Particles%wpv(i)%is_mapped) &
                      nb_wpv = nb_wpv + 1
              ENDDO
              ALLOCATE(wpv_l(nb_wpv),STAT=info)
              nb_wpv = 0
              DO i=1,Particles%max_wpvid
                 IF (Particles%wpv(i)%is_mapped) THEN
                    nb_wpv = nb_wpv + 1
                    wpv_l(nb_wpv) = i
                 ENDIF
              ENDDO
           ENDIF

           ! vector field property
           IF(PRESENT(wpv_field_list)) THEN
              nb_wpv_field=SIZE(wpv_field_list)
              ALLOCATE(wpv_field_l(nb_wpv_field),STAT=info)
              wpv_field_l=wpv_field_list
              DO i=1,nb_wpv_field
                 IF (wpv_field_l(i).GT.Particles%max_wpvid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'vector field property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (.NOT.Particles%wpv(wpv_field_l(i))%is_mapped) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              ! dont print anything
              nb_wpv_field = 0
           ENDIF


#ifdef __MPI
           ! write parallel file
           IF (ppm_rank .EQ. 0) THEN
              WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.pvtp'
              OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
                   IOSTAT=info, ACTION='WRITE')
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 WRITE(errtxt,'(2A)') 'Failed to open file: ', &
                      scratch(1:LEN_TRIM(scratch))
                 CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
                 GOTO 9999
              END IF
#define VTK_FILE_TYPE "PPolyData"
#define VTK_PARALLEL
#include "vtk/print_header.f"
              WRITE(iUnit,'(A)') "    <PPointData>"
              IF (nvlist) THEN
              WRITE(iUnit,'(3A)') "      <PDataArray Name='nvlist' type='Float64' />"
              END IF
              DO i=1,nb_wpi
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   Particles%wpi(wpi_l(i))%name &
                   (1:LEN_TRIM(Particles%wpi(wpi_l(i))%name)), &
                   "' type='Float64' />"
              END DO
              DO i=1,nb_wps
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   Particles%wps(wps_l(i))%name &
                   (1:LEN_TRIM(Particles%wps(wps_l(i))%name)), &
                   "' type='Float64' />"
              END DO
              DO i=1,nb_wpv
                 nd = Particles%wpv(i)%lda
                 DO l=1,nd
                    WRITE(scratch,'(A,A,I0)') TRIM(Particles%wpv(wpv_l(i))%name), '_', l
                    WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                        scratch(1:LEN_TRIM(scratch)), "' type='Float64' />"
                 END DO
              END DO
              DO i=1,nb_wpv_field
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   Particles%wpv(wpv_field_l(i))%name &
                   (1:LEN_TRIM(Particles%wpv(wpv_field_l(i))%name)), &
                   "' type='Float64' />"
              END DO
              WRITE(iUnit,'(A)') "    </PPointData>"              
              WRITE(iUnit,'(A)') "    <PPoints>"
              WRITE(iUnit,'(A)') "      <PDataArray NumberOfComponents='3' type='Float64' />"
              WRITE(iUnit,'(A)') "    </PPoints>"
              ! find the basename of the file
              DO i=0,ppm_nproc-1
                 WRITE(iUnit,'(A,A,A,I0,A)') "    <Piece Source='",     &
                      fname(INDEX(fname, '/', .true.)+1:LEN_TRIM(fname)), &
                      ".", i, ".vtp' />"
              END DO
              ! close
#include "vtk/print_end_header.f"
              CLOSE(iUnit)
           END IF
           ! append rank to name
           WRITE(scratch,'(A,A,I0,A)') fname(1:LEN_TRIM(fname)), &
                                       '.', ppm_rank, '.vtp'
#else
           WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.vtp'
#endif

           ! open output file
           OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
                IOSTAT=info, ACTION='WRITE')
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              WRITE(errtxt,'(2A)') 'Failed to open file: ', &
                    scratch(1:LEN_TRIM(scratch))
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
           END IF

           ! write data
           IF (ghosts) THEN
              N = Particles%Mpart
           ELSE
              N = Particles%Npart
           END IF

           ! write header
#define VTK_FILE_TYPE "PolyData"
#define VTK_NPOINTS N
#define VTK_NVERTS  N
#include "vtk/print_header.f"

           ! print properties
           IF (nvlist .OR. nb_wpi .GT. 0 .OR. nb_wps .GT. 0 .OR. nb_wpv .GT. 0) THEN

              ! print names
              WRITE(iUnit,'(A)',advance='no') "      <PointData" 
              IF (nvlist .OR. nb_wpi .GT. 0) THEN
                 WRITE(iUnit,'(A)',advance='no') " Integers='"
              ENDIF
              IF (nvlist) THEN
                 WRITE(iUnit,'(A)',advance='no') "nvlist"
                 IF (nb_wpi .GT. 0) WRITE(iUnit,'(A)',advance='no') " "
              END IF
              IF (nb_wpi .GT. 0) THEN
                 DO i=1,nb_wpi
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wpi(wpi_l(i))%name &
                         (1:LEN_TRIM(Particles%wpi(wpi_l(i))%name))
                    IF (i .LT. nb_wpi) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              IF (nb_wps .GT. 0) THEN
                 IF (nvlist .OR. nb_wpi .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') "'"
                 WRITE(iUnit,'(A)',advance='no') " Scalars='"
                 DO i=1,nb_wps
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wps(wps_l(i))%name &
                         (1:LEN_TRIM(Particles%wps(wps_l(i))%name))
                    IF (i .LT. nb_wps) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              IF (nb_wpv .GT. 0) THEN
                 IF (nvlist .OR. nb_wpi .GT. 0 .OR. nb_wps .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') "'"
                 WRITE(iUnit,'(A)',advance='no') " Vectors='"
                 DO i=1,nb_wpv
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wpv(wpv_l(i))%name &
                         (1:LEN_TRIM(Particles%wpv(wpv_l(i))%name))
                    IF (i .LT. nb_wpv) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
                 IF (nb_wpv .GT. 0 .AND. nb_wpv_field .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') " "
                 DO i=1,nb_wpv_field
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wpv(wpv_field_l(i))%name &
                         (1:LEN_TRIM(Particles%wpv(wpv_field_l(i))%name))
                    IF (i .LT. nb_wpv_field) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              WRITE(iUnit,'(A)') "'>"

              ! property values
              IF (nvlist) THEN
                 wpi => Particles%nvlist
#define VTK_NAME "nvlist"
#define VTK_TYPE "Float64"
#define VTK_INTEGER wpi
#include "vtk/print_data_array.f"
                 wpi => NULL()
              END IF
              DO k=1,nb_wpi
                 wpi => get_wpi(Particles,wpi_l(k),with_ghosts=ghosts)
#define VTK_NAME Particles%wpi(wpi_l(k))%name
#define VTK_TYPE "Float64"
#define VTK_INTEGER wpi
#include "vtk/print_data_array.f"
                 wpi => set_wpi(Particles,wpi_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wps
                 wp => get_wps(Particles,wps_l(k),with_ghosts=ghosts)
#define VTK_NAME Particles%wps(wps_l(k))%name
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
                 wp => set_wps(Particles,wps_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wpv
                 xp => get_wpv(Particles,wpv_l(k),with_ghosts=ghosts)
                 nd = SIZE(xp,1)
                 DO l=1,nd
                    WRITE(scratch,'(A,A,I0)') TRIM(Particles%wpv(wpv_l(k))%name), '_', l
                    wp => xp(l,:)
#define VTK_NAME scratch
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
                    wp => NULL()
                 END DO
                 xp => set_wpv(Particles,wpv_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wpv_field
                 xp => get_wpv(Particles,wpv_field_l(k),with_ghosts=ghosts)
                 nd = SIZE(xp,1)
#define VTK_NAME Particles%wpv(wpv_field_l(k))%name
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
                 xp => set_wpv(Particles,wpv_field_l(k),read_only=.TRUE.)
              END DO
              WRITE(iUnit,'(A)') "      </PointData>"
           END IF

           ! print point coordinates
           WRITE(iUnit,'(A)') "      <Points>"
           xp => get_xp(Particles,with_ghosts=ghosts)
           nd = SIZE(xp,1)
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
           xp => set_xp(Particles,read_only=.TRUE.)
           WRITE(iUnit,'(A)') "      </Points>"

           ! create a vertex for every point
           WRITE(iUnit,'(A)') "      <Verts>"
           ! connectivity
           N = N - 1
#define VTK_RANGE N
#define VTK_RANGE_START 0
#define VTK_NAME "connectivity"
#define VTK_TYPE "Int32"
#include "vtk/print_data_array.f"
           ! offsets
           N = N + 1
#define VTK_RANGE N
#define VTK_NAME "offsets"
#define VTK_TYPE "Int32"
#include "vtk/print_data_array.f"
           WRITE(iUnit,'(A)') "      </Verts>"

           ! close
#include "vtk/print_end_header.f"
           ! close file
9999       CONTINUE
           CLOSE(iUnit)
           CALL substop(caller,t0,info)
         END SUBROUTINE DTYPE(ppm_vtk_particle_cloud)

#undef __KIND
#undef DTYPE
#undef DEFINE_MK
