         SUBROUTINE DTYPE(ppm_vtk_particles)(filename, Pc, info, &
              step, with_ghosts, with_nvlist, Fields)
           !--------------------------------------------------------------------
           !  Arguments
           !--------------------------------------------------------------------
           DEFINE_MK()
           CHARACTER(LEN=*),                INTENT(IN   ) :: filename
           CLASS(DTYPE(ppm_t_particles)_),  INTENT(INOUT) :: Pc
           INTEGER,                         INTENT(  OUT) :: info
           INTEGER,               OPTIONAL, INTENT(IN   ) :: step
           LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_ghosts
           LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_nvlist
           CLASS(ppm_v_main_abstr),OPTIONAL               :: Fields
           !!!list of fields to print out (default = ALL)

           !INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wp_list
           !list of properties to print out (default = ALL)
           !INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wp_field_list
           !list of properties to print out as vector field (default = NONE)
           !--------------------------------------------------------------------
           !  Variables
           !--------------------------------------------------------------------
           CHARACTER(LEN=ppm_char)              :: scratch
           CHARACTER(LEN=ppm_char)              :: fname
           INTEGER                              :: i,j,k,l,nd,N,ii
           INTEGER                              :: nb_wps, nb_wpv, nb_wp_field
           INTEGER                              :: nb_wp,nb_wpi
           INTEGER, DIMENSION(:),   ALLOCATABLE :: wpi_l, wps_l, wpv_l, wp_field_l
           LOGICAL                              :: ghosts
           LOGICAL                              :: nvlist
           REAL(MK), DIMENSION(:,:),POINTER     :: xp  => NULL()
           REAL(MK), DIMENSION(:),  POINTER     :: wp  => NULL()
           REAL(MK), DIMENSION(:,:),  POINTER   :: wp2d  => NULL()
           INTEGER, DIMENSION(:),   POINTER     :: wpi  => NULL()

           CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()
           CLASS(ppm_t_discr_data),       POINTER :: discr_data => NULL()
           CLASS(ppm_t_field_),           POINTER :: field => NULL()
           CLASS(ppm_t_main_abstr),       POINTER :: el => NULL()

           start_subroutine("ppm_vtk_particles")
           
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

           nb_wpi=0
           nb_wps=0
           nb_wpv=0

           ! create the list of properties to print
           ! integer property
           !IF(PRESENT(wp_list)) THEN
              !nb_wp=SIZE(wp_list)

              !ALLOCATE(wpi_l(nb_wp),wps_l(nb_wp),wpv_l(nb_wp),STAT=info)
              !!find how many integer/scalar/vector properties are to be
              !!printed out
              !DO i=1,nb_wp
                 !ii = wp_list(i)
                 !IF (.NOT.Pc%props%exists(ii)) THEN
                    !info = ppm_error_error
                    !CALL ppm_error(ppm_err_argument,caller,   &
                 !&  'invalid property for printout. Data may not be allocated',&
                         !&  __LINE__,info)
                    !GOTO 9999
                 !ENDIF
                 !prop => Pc%props%vec(ii)%t
                 !IF (.NOT.prop%flags(ppm_ppt_partial)) THEN
                    !info = ppm_error_error
                    !CALL ppm_error(ppm_err_argument,caller,   &
                         !&  'trying to printout a property that is not mapped &
                         !& to the particles',&
                         !&  __LINE__,info)
                    !GOTO 9999
                 !ENDIF
                 !SELECT CASE (prop%data_type)
                 !CASE (ppm_type_int)
                     !IF (prop%lda.NE.1) THEN
                         !info = ppm_error_error
                         !CALL ppm_error(ppm_err_argument,caller,   &
                     !&  'no support for printout of vector integer properties',&
                             !&  __LINE__,info)
                         !GOTO 9999
                     !ENDIF
                     !nb_wpi = nb_wpi + 1
                     !wpi_l(nb_wpi) = ii

                 !CASE (ppm_type_real)
                     !IF (prop%lda.EQ.1) THEN
                         !nb_wps = nb_wps + 1
                         !wps_l(nb_wps) = ii
                     !ELSE
                         !nb_wpv = nb_wpv + 1
                         !wpv_l(nb_wpv) = ii
                     !ENDIF
                 !CASE DEFAULT
                     !info = ppm_error_error
                     !CALL ppm_error(ppm_err_argument,caller,   &
                 !&  'no support for vtk printout of data other than real or int',&
                         !&  __LINE__,info)
                     !GOTO 9999
                 !END SELECT
                 !prop => NULL()
              !ENDDO

           !ELSE
              !nb_wp = 0
           !ENDIF


           !! vector field property
           !IF(PRESENT(wp_field_list)) THEN
              !nb_wp_field=SIZE(wp_field_list)
              !ALLOCATE(wp_field_l(nb_wp_field),STAT=info)
              !wp_field_l=wp_field_list

              !!Check that the listed properties are indeed vector fields,
              !! are allocated and mapped to particles
              !DO i=1,nb_wp_field
                 !ii = wp_field_list(i)
                 !IF (.NOT.Pc%props%exists(ii)) THEN
                    !info = ppm_error_error
                    !CALL ppm_error(ppm_err_argument,caller,   &
                 !&  'invalid property for printout. Data may not be allocated',&
                         !&  __LINE__,info)
                    !GOTO 9999
                 !ENDIF
                 !prop => Pc%props%vec(ii)%t
                 !IF (.NOT.prop%flags(ppm_ppt_partial)) THEN
                    !info = ppm_error_error
                    !CALL ppm_error(ppm_err_argument,caller,   &
                         !&  'trying to printout a property that is not mapped &
                         !& to the particles',&
                         !&  __LINE__,info)
                    !GOTO 9999
                 !ENDIF
                 !IF (.NOT. prop%data_type.EQ.ppm_type_real .OR. &
                        !prop%lda.LE.1) THEN 
                         !info = ppm_error_error
                         !CALL ppm_error(ppm_err_argument,caller,   &
                     !&  'property should be a real vector field',&
                             !&  __LINE__,info)
                         !GOTO 9999
                 !ENDIF
                 !prop => NULL()
              !ENDDO
           !ELSE
              !! dont print anything
              !nb_wp_field = 0
           !ENDIF

           !IF (.NOT. PRESENT(wp_list) .AND. .NOT. PRESENT(wp_field_list)) THEN
               !!printout all properties i that are mapped
               !nb_wp = Pc%props%nb
               !ALLOCATE(wpi_l(nb_wp),wps_l(nb_wp),wpv_l(nb_wp),STAT=info)
               !prop => Pc%props%begin()
               !DO WHILE (ASSOCIATED(prop))
                   !IF (prop%flags(ppm_ppt_partial)) THEN
                       !SELECT CASE (prop%data_type)
                       !CASE (ppm_type_int)
                           !IF (prop%lda.EQ.1) THEN
                               !nb_wpi = nb_wpi + 1
                               !wpi_l(nb_wpi) = Pc%props%iter_id
                           !ENDIF

                       !CASE (ppm_type_real)
                           !IF (prop%lda.EQ.1) THEN
                               !nb_wps = nb_wps + 1
                               !wps_l(nb_wps) = Pc%props%iter_id
                           !ELSE
                               !nb_wpv = nb_wpv + 1
                               !wpv_l(nb_wpv) = Pc%props%iter_id
                           !ENDIF
                       !CASE DEFAULT
                           !!not a supported type for printout (yet)
                       !END SELECT
                       !prop => Pc%props%next()
                   !ENDIF
               !ENDDO
           !ENDIF

           IF (PRESENT(Fields)) THEN
               el => Fields%begin()
               DO WHILE (ASSOCIATED(field))
                   SELECT TYPE(field => el)
                   CLASS IS(ppm_t_field_)
                       !hack, so that Pc%props%iter_id is now the id
                       !of the discretization of field in Pc.
                       CALL field%get_discr(Pc,discr_data,info)
                       or_fail("could not get discr data for this field")
                       SELECT CASE(field%lda)
                       CASE (ppm_type_int)
                           nb_wpi = nb_wpi + 1
                           wpi_l(nb_wpi) = Pc%props%iter_id
                       CASE (ppm_type_real)
                           IF (field%lda.EQ.1) THEN
                               nb_wps = nb_wps + 1
                               wps_l(nb_wps) = Pc%props%iter_id
                           ELSE
                               nb_wpv = nb_wpv + 1
                               wpv_l(nb_wpv) = Pc%props%iter_id
                           ENDIF
                       CASE DEFAULT
                               !not a supported type for printout (yet)
                       END SELECT
                   END SELECT
                   el => Fields%next()
               ENDDO
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
                  prop => Pc%props%vec(wpi_l(i))%t
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                      prop%name (1:LEN_TRIM(prop%name)), &
                      "' type='Float64' />"
              END DO
              DO i=1,nb_wps
                  prop => Pc%props%vec(wps_l(i))%t
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                      prop%name (1:LEN_TRIM(prop%name)), &
                      "' type='Float64' />"
              END DO
              DO i=1,nb_wpv
                  prop => Pc%props%vec(wpv_l(i))%t
                  nd = prop%lda
                  DO l=1,nd
                      WRITE(scratch,'(A,A,I0)') TRIM(prop%name), '_', l
                      WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                          scratch(1:LEN_TRIM(scratch)), "' type='Float64' />"
                  END DO
              END DO
              DO i=1,nb_wp_field
                  prop => Pc%props%vec(wp_field_l(i))%t
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   prop%name (1:LEN_TRIM(prop%name)), &
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
              N = Pc%Mpart
           ELSE
              N = Pc%Npart
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
                    prop => Pc%props%vec(wpi_l(i))%t
                    WRITE(iUnit,'(A)',advance='no') &
                         prop%name (1:LEN_TRIM(prop%name))
                    IF (i .LT. nb_wpi) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              IF (nb_wps .GT. 0) THEN
                 IF (nvlist .OR. nb_wpi .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') "'"
                 WRITE(iUnit,'(A)',advance='no') " Scalars='"
                 DO i=1,nb_wps
                    prop => Pc%props%vec(wps_l(i))%t
                    WRITE(iUnit,'(A)',advance='no') &
                         prop%name (1:LEN_TRIM(prop%name))
                    IF (i .LT. nb_wps) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              IF (nb_wpv .GT. 0) THEN
                 IF (nvlist .OR. nb_wpi .GT. 0 .OR. nb_wps .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') "'"
                 WRITE(iUnit,'(A)',advance='no') " Vectors='"
                 DO i=1,nb_wpv
                    prop => Pc%props%vec(wpv_l(i))%t
                    WRITE(iUnit,'(A)',advance='no') &
                         prop%name (1:LEN_TRIM(prop%name))
                    IF (i .LT. nb_wpv) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
                 IF (nb_wpv .GT. 0 .AND. nb_wp_field .GT. 0) &
                      WRITE(iUnit,'(A)',advance='no') " "
                 DO i=1,nb_wp_field
                    prop => Pc%props%vec(wp_field_l(i))%t
                    WRITE(iUnit,'(A)',advance='no') &
                         prop%name (1:LEN_TRIM(prop%name))
                    IF (i .LT. nb_wp_field) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              WRITE(iUnit,'(A)') "'>"

              ! property values
              IF (nvlist) THEN
                 CALL Pc%get_nvlist(nvlist=wpi,info=info)
                    or_fail("could not access neighbour list")
#define VTK_NAME "nvlist"
#define VTK_TYPE "Float64"
#define VTK_INTEGER wpi
#include "vtk/print_data_array.f"
                 wpi => null()
              end if
              do k=1,nb_wpi
                 prop => pc%props%vec(wpi_l(k))%t
                 CALL Pc%get(wpi,wpi_l(k),with_ghosts=ghosts)
#define VTK_NAME prop%name
#define VTK_TYPE "Float64"
#define VTK_INTEGER wpi
#include "vtk/print_data_array.f"
                 CALL Pc%set(wpi,wpi_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wps
                 prop => Pc%props%vec(wps_l(k))%t
                 CALL Pc%get(wp,wps_l(k),with_ghosts=ghosts)
#define VTK_NAME prop%name
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
                 CALL Pc%set(wp,wps_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wpv
                 prop => Pc%props%vec(wpv_l(k))%t
                 CALL Pc%get(wp2d,wpv_l(k),with_ghosts=ghosts)
                 nd = SIZE(wp2d,1)
                 DO l=1,nd
                    WRITE(scratch,'(A,A,I0)') TRIM(prop%name), '_', l
                    wp => wp2d(l,:)
#define VTK_NAME scratch
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
                    wp => NULL()
                 END DO
                 CALL Pc%set(wp2d,wpv_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wp_field
                 prop => Pc%props%vec(wp_field_l(k))%t
                 CALL Pc%get(wp2d,wp_field_l(k),with_ghosts=ghosts)
                 nd = SIZE(wp2d,1)
#define VTK_NAME prop%name
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR wp2d
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
                 CALL Pc%set(wp2d,wp_field_l(k),read_only=.TRUE.)
              END DO
              WRITE(iUnit,'(A)') "      </PointData>"
           END IF

           ! print point coordinates
           WRITE(iUnit,'(A)') "      <Points>"
           CALL Pc%get_xp(xp,info,with_ghosts=ghosts)
                or_fail("get_xp")
           nd = SIZE(xp,1)
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
           CALL Pc%set_xp(xp,info,read_only=.TRUE.)
                or_fail("set_xp")
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
           end_subroutine()
           CLOSE(iUnit)
         END SUBROUTINE DTYPE(ppm_vtk_particles)

#undef __KIND
#undef DTYPE
#undef DEFINE_MK
